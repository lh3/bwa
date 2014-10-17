#include <sys/types.h>
#include <sys/mman.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include "bwa.h"

int bwa_shm_stage(bwaidx_t *idx, const char *hint, const char *_tmpfn)
{
	const char *name;
	uint8_t *shm, *shm_idx;
	uint16_t *cnt;
	int shmid, to_init = 0, l;
	char path[PATH_MAX + 1], *tmpfn = (char*)_tmpfn;

	if (hint == 0 || hint[0] == 0) return -1;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;

	if ((shmid = shm_open("/bwactl", O_RDWR, 0)) < 0) {
		shmid = shm_open("/bwactl", O_CREAT|O_RDWR|O_EXCL, 0644);
		to_init = 1;
	}
	if (shmid < 0) return -1;
	ftruncate(shmid, BWA_CTL_SIZE);
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	if (to_init) {
		memset(shm, 0, BWA_CTL_SIZE);
		cnt[1] = 4;
	}

	if (idx->mem == 0) bwa_idx2mem(idx);

	if (tmpfn) {
		FILE *fp;
		if ((fp = fopen(tmpfn, "wb")) != 0) {
			int64_t rest = idx->l_mem;
			while (rest > 0) {
				int64_t l = rest < 0x1000000? rest : 0x1000000;
				rest -= fwrite(&idx->mem[idx->l_mem - rest], 1, l, fp);
			}
			fclose(fp);
			free(idx->mem); idx->mem = 0;
		} else {
			fprintf(stderr, "[W::%s] fail to create the temporary file. Option '-f' is ignored.\n", __func__);
			tmpfn = 0;
		}
	}

	strcat(strcpy(path, "/bwaidx-"), name);
	if ((shmid = shm_open(path, O_CREAT|O_RDWR|O_EXCL, 0644)) < 0) {
		shm_unlink(path);
		perror("shm_open()");
		return -1;
	}
	l = 8 + strlen(name) + 1;
	if (cnt[1] + l > BWA_CTL_SIZE) return -1;
	memcpy(shm + cnt[1], &idx->l_mem, 8);
	memcpy(shm + cnt[1] + 8, name, l - 8);
	cnt[1] += l; ++cnt[0];
	ftruncate(shmid, idx->l_mem);
	shm_idx = mmap(0, idx->l_mem, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	if (tmpfn) {
		FILE *fp;
		fp = fopen(tmpfn, "rb");
		int64_t rest = idx->l_mem;
		while (rest > 0) {
			int64_t l = rest < 0x1000000? rest : 0x1000000;
			rest -= fread(&shm_idx[idx->l_mem - rest], 1, l, fp);
		}
		fclose(fp);
		unlink(tmpfn);
	} else {
		memcpy(shm_idx, idx->mem, idx->l_mem);
		free(idx->mem);
	}
	bwa_mem2idx(idx->l_mem, shm_idx, idx);
	idx->is_shm = 1;
	return 0;
}

bwaidx_t *bwa_idx_load_from_shm(const char *hint)
{
	const char *name;
	uint8_t *shm, *shm_idx;
	uint16_t *cnt, i;
	char *p, path[PATH_MAX + 1];
	int shmid;
	int64_t l_mem;
	bwaidx_t *idx;

	if (hint == 0 || hint[0] == 0) return 0;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;
	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return 0;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	if (cnt[0] == 0) return 0;
	for (i = 0, p = (char*)(shm + 4); i < cnt[0]; ++i) {
		memcpy(&l_mem, p, 8); p += 8;
		if (strcmp(p, name) == 0) break;
		p += strlen(p) + 1;
	}
	if (i == cnt[0]) return 0;

	strcat(strcpy(path, "/bwaidx-"), name);
	if ((shmid = shm_open(path, O_RDONLY, 0)) < 0) return 0;
	shm_idx = mmap(0, l_mem, PROT_READ, MAP_SHARED, shmid, 0);
	idx = calloc(1, sizeof(bwaidx_t));
	bwa_mem2idx(l_mem, shm_idx, idx);
	idx->is_shm = 1;
	return idx;
}

int bwa_shm_test(const char *hint)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	const char *name;

	if (hint == 0 || hint[0] == 0) return 0;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;
	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return 0;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		if (strcmp(p + 8, name) == 0) return 1;
		p += strlen(p) + 9;
	}
	return 0;
}

int bwa_shm_list(void)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		int64_t l_mem;
		memcpy(&l_mem, p, 8); p += 8;
		printf("%s\t%ld\n", p, (long)l_mem);
		p += strlen(p) + 1;
	}
	return 0;
}

int bwa_shm_destroy(void)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	char path[PATH_MAX + 1];

	if ((shmid = shm_open("/bwactl", O_RDONLY, 0)) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		int64_t l_mem;
		memcpy(&l_mem, p, 8); p += 8;
		strcat(strcpy(path, "/bwaidx-"), p);
		shm_unlink(path);
		p += strlen(p) + 1;
	}
	munmap(shm, BWA_CTL_SIZE);
	shm_unlink("/bwactl");
	return 0;
}

int main_shm(int argc, char *argv[])
{
	int c, to_list = 0, to_drop = 0, ret = 0;
	char *tmpfn = 0;
	while ((c = getopt(argc, argv, "ldf:")) >= 0) {
		if (c == 'l') to_list = 1;
		else if (c == 'd') to_drop = 1;
		else if (c == 'f') tmpfn = optarg;
	}
	if (optind == argc && !to_list && !to_drop) {
		fprintf(stderr, "\nUsage: bwa shm [-d|-l] [-f tmpFile] [idxbase]\n\n");
		fprintf(stderr, "Options: -d       destroy all indices in shared memory\n");
		fprintf(stderr, "         -l       list names of indices in shared memory\n");
		fprintf(stderr, "         -f FILE  temporary file to reduce peak memory\n\n");
		return 1;
	}
	if (optind < argc && (to_list || to_drop)) {
		fprintf(stderr, "[E::%s] open -l or -d cannot be used when 'idxbase' is present\n", __func__);
		return 1;
	}
	if (optind < argc) {
		if (bwa_shm_test(argv[optind]) == 0) {
			bwaidx_t *idx;
			idx = bwa_idx_load_from_disk(argv[optind], BWA_IDX_ALL);
			if (bwa_shm_stage(idx, argv[optind], tmpfn) < 0) {
				fprintf(stderr, "[E::%s] failed to stage the index in shared memory\n", __func__);
				ret = 1;
			}
			bwa_idx_destroy(idx);
		} else fprintf(stderr, "[M::%s] index '%s' is already in shared memory\n", __func__, argv[optind]);
	}
	if (to_list) bwa_shm_list();
	if (to_drop) bwa_shm_destroy();
	return ret;
}
