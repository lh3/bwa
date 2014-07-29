#pragma once

#include <stdint.h>
#include <cstdlib>
#include <vector>
#include <iostream>

#include "tbb/pipeline.h"
#include "tbb/tick_count.h"
#include "tbb/pipeline.h"

#include "bwa.h"
#include "kseq.h"

extern "C"
{
	KSEQ_DECLARE(gzFile)
}
extern char *bwa_pg;

extern unsigned char nst_nt4_table[256];

extern "C" void *kopen(const char *fn, int *_fd);
extern "C" int kclose(void *a);
extern "C" int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);
extern "C" mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq);
extern "C" void mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id); // IMPORTANT: must run mem_sort_and_dedup() before calling this function
extern "C" void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);


namespace bwa
{
	class BSeq
	{
	public:
		BSeq(bseq1_t* seqs, uint32_t size, uint64_t processed): seqs_(seqs), size_(size), processed_(processed)
		{
		}

		~BSeq()
		{
			for (size_t i = 0; i < size_; i++)
			{
				free(seqs_[i].name);
				free(seqs_[i].comment);
				free(seqs_[i].seq);
				free(seqs_[i].qual);
				free(seqs_[i].sam);
			}
			free(seqs_);
			size_=0;
		}

		bseq1_t* data()
		{
			return seqs_;
		}

		uint32_t size()
		{
			return size_;
		}

		uint64_t processed()
		{
			return processed_;
		}

	private:
		bseq1_t* seqs_;
		size_t size_;
		const uint64_t processed_;
	}; 

	template<typename FILTER, typename IN_TYPE, typename OUT_TYPE>
	class LogFilter
	{
	public:
		LogFilter(const std::string& name, const FILTER& filter): name_(name), filter_(filter)
		{
			std::clog<<name_<<std::endl;
		}

		OUT_TYPE operator()(IN_TYPE in) const
		{
			tbb::tick_count t=tbb::tick_count::now();
			std::clog<<name_<<"()"<<std::endl;
			OUT_TYPE out=filter_(in);
			std::clog<<name_<<"() time="<<(tbb::tick_count::now()-t).seconds()<<std::endl;
			return out;
		}

		const FILTER* operator->() const
		{
			return &filter_;
		}

		~LogFilter()
		{
			std::clog<<"~"<<name_<<std::endl;
		}

	private:
		const std::string name_;
		const FILTER filter_;
	};

	template<typename FILTER, typename IN_TYPE>
	class LogFilter<FILTER, IN_TYPE, void>
	{
	public:
		LogFilter(const std::string& name, const FILTER& filter): name_(name), filter_(filter)
		{
			std::clog<<name_<<std::endl;
		}

		void operator()(IN_TYPE in) const
		{
			tbb::tick_count t=tbb::tick_count::now();
			std::clog<<name_<<"()"<<std::endl;
			filter_(in);
			std::clog<<name_<<"() time="<<(tbb::tick_count::now()-t).seconds()<<std::endl;
		}

		const FILTER* operator->() const
		{
			return &filter_;
		}

		~LogFilter()
		{
			std::clog<<"~"<<name_<<std::endl;
		}

	private:
		const std::string name_;
		const FILTER filter_;
	};

	class ReadBSeqFilter
	{
	public:
		ReadBSeqFilter(std::vector<kseq_t*> ks, const mem_opt_t& opt, bool copy_comment): ks_(ks), opt_(opt), copy_comment_(copy_comment)
		{
			processed_=0;
			threads_=0;
		}

		~ReadBSeqFilter()
		{
		}

		BSeq* operator()(tbb::flow_control& fc) const
		{
			int n_seqs;
			bseq1_t* seqs;
			if(seqs=bseq_read(opt_.chunk_size,&n_seqs,ks_[0],ks_[1]))
			{
				int64_t size = 0;
				if ((opt_.flag & MEM_F_PE) && (n_seqs&1) == 1) {
					if (bwa_verbose >= 2)
						fprintf(stderr, "[W::%s] odd number of reads in the PE mode; last read dropped\n", __func__);
					n_seqs = n_seqs>>1<<1;
				}
				if (!copy_comment_)
					for (int i = 0; i < n_seqs; ++i) {
						free(seqs[i].comment); seqs[i].comment = 0;
					}
					for (int i = 0; i < n_seqs; ++i) size += seqs[i].l_seq;
					if (bwa_verbose >= 3)
						fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, n_seqs, (long)size);
			}
			else
			{
				fc.stop();
				return nullptr;
			}

			uint64_t processed=processed_;
			processed_+=n_seqs;
			return new BSeq(seqs,n_seqs,processed);
		}

	private:
		const mem_opt_t& opt_;
		std::vector<kseq_t*> ks_;
		const bool copy_comment_;
		mutable uint64_t processed_;
		mutable uint32_t threads_;
	};

	class WriteBSeqFilter
	{
	public:
		WriteBSeqFilter(std::ostream& out, const bntseq_t& bns, const std::string& rg_line): out_(out), bns_(bns), rg_line_(rg_line)
		{
		}

		~WriteBSeqFilter()
		{
		}

		void write_hdr() const
		{
			for (int i = 0; i < bns_.n_seqs; ++i)
			{
				//err_printf("@SQ\tSN:%s\tLN:%d\n", bns_->anns[i].name, bns_->anns[i].len);
				out_<<"@SQ\tSN:"<<bns_.anns[i].name<<"\tLN:"<<bns_.anns[i].len<<"\n";
			}
			if (!rg_line_.empty())
			{
				//err_printf("%s\n", rg_line_);
				out_<<rg_line_<<"\n";
			}
			out_<<bwa_pg<<"\n";
		}

		void operator()(BSeq* seqs) const
		{
			for(size_t i=0; i<seqs->size(); i++)
			{
				out_<<(seqs->data())[i].sam;
			}
			delete seqs;
		}

	private:
		const bntseq_t& bns_;
		const std::string& rg_line_;
		std::ostream& out_;
	};


	class MemParallelFilter
	{
	public:
		//const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0)
		MemParallelFilter(const bwt_t& bwt,  const bntseq_t& bns, const uint8_t* pacseq, const mem_opt_t& opt, const mem_pestat_t* pes0):
			bwt_(bwt), bns_(bns), pacseq_(pacseq), opt_(opt), pes0_(pes0)
		{
		}

		~MemParallelFilter()
		{

		}

		BSeq* operator()(BSeq* seqs) const
		{
			mem_alnreg_v* regs=reinterpret_cast<mem_alnreg_v*>(malloc(seqs->size() * sizeof(mem_alnreg_v)));
			mem_pestat_t pes[4];

			int n_seqs=(opt_.flag&MEM_F_PE)? seqs->size()>>1 : seqs->size();


			for(size_t i=0; i<n_seqs; i++)
			{
				//worker1
				if(!(opt_.flag&MEM_F_PE))
				{
					if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", seqs->data()[i].name);
					regs[i] = mem_align1_core(&opt_, &bwt_, &bns_, pacseq_, seqs->data()[i].l_seq, seqs->data()[i].seq);
				}
				else
				{
					if(bwa_verbose >= 4) printf("=====> Processing read '%s'/1 <=====\n", seqs->data()[i<<1|0].name);
					regs[i<<1|0] = mem_align1_core(&opt_, &bwt_, &bns_, pacseq_, seqs->data()[i<<1|0].l_seq, seqs->data()[i<<1|0].seq);
					if(bwa_verbose >= 4) printf("=====> Processing read '%s'/2 <=====\n", seqs->data()[i<<1|1].name);
					regs[i<<1|1] = mem_align1_core(&opt_, &bwt_, &bns_, pacseq_, seqs->data()[i<<1|1].l_seq, seqs->data()[i<<1|1].seq);
				}
			}

			//TODO: behaviour should be the same as original version with one thread, in other case this filter should be splitted to three...
			if(opt_.flag&MEM_F_PE)
			{ // infer insert sizes if not provided
				if(pes0_) memcpy(pes, pes0_, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
				else mem_pestat(&opt_, bns_.l_pac, seqs->size(), regs, pes); // otherwise, infer the insert size distribution from data
			}

			for(size_t i=0; i<n_seqs; i++)
			{
				//worker 2
				if(!(opt_.flag&MEM_F_PE))
				{
					if(bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", seqs->data()[i].name);
					mem_mark_primary_se(&opt_, regs[i].n, regs[i].a, /*w->n_processed +*/ seqs->processed()+i);
					mem_reg2sam_se(&opt_, &bns_, pacseq_, &seqs->data()[i], &regs[i], 0, 0);
					free(regs[i].a);
				}
				else
				{
					if(bwa_verbose >= 4) printf("=====> Finalizing read pair '%s' <=====\n", seqs->data()[i<<1|0].name);
					mem_sam_pe(&opt_, &bns_, pacseq_, &pes[0], /*(w->n_processed>>1) +*/(seqs->processed()>>1)+i, &seqs->data()[i<<1], &regs[i<<1]);
					free(regs[i<<1|0].a); free(regs[i<<1|1].a);
				}
			}
			free(regs);
			return seqs;
		}

	private:
		const bntseq_t& bns_;
		const mem_pestat_t* pes0_;
		const mem_opt_t& opt_;
		const bwt_t& bwt_;
		const uint8_t* pacseq_;
	};
}
