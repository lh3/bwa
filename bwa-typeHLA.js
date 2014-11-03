var fuzzy = 3, min_qal = 100, drop_thres = 7;

var file = arguments.length == 0? new File () : new File(arguments[0]);
var buf = new Bytes();
var re_cigar = /(\d+)([MIDSH])/g;

var len = {}, list = [], gcnt = [];
while (file.readline(buf) >= 0) {
	var m, mm, line = buf.toString();
	var t = line.split("\t");
	if (t[0].charAt(0) == '@') {
		if (t[0] == '@SQ' && (m = /LN:(\d+)/.exec(line)) != null && (mm = /SN:(\S+)/.exec(line)) != null)
			len[mm[1]] = parseInt(m[1]);
		continue;
	}
	var gene = null, exon = null;
	if ((m = /^(HLA-[^\s_]+)_(\d+)/.exec(t[0])) != null) {
		gene = m[1], exon = parseInt(m[2]) - 1;
		if (gcnt[exon] == null) gcnt[exon] = {};
		gcnt[exon][gene] = true;
	}
	if (gene == null || exon == null || t[2] == '*') continue;
	var x = 0, ts = parseInt(t[3]) - 1, te = ts, clip = [0, 0];
	while ((m = re_cigar.exec(t[5])) != null) {
		var l = parseInt(m[1]);
		if (m[2] == 'M') x += l, te += l;
		else if (m[2] == 'I') x += l;
		else if (m[2] == 'D') te += l;
		else if (m[2] == 'S' || m[2] == 'H') clip[x==0?0:1] = l;
	}
	if (x < min_qal && clip[0] + clip[1] > 0) continue;
	var tl = len[t[2]];
	var left  = ts < clip[0]? ts : clip[0];
	var right = tl - te < clip[1]? tl - te : clip[1];
	var nm = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : 0;
	list.push([t[2], gene, exon, ts, te, nm, left + right]);
}

buf.destroy();
file.close();

// identify the primary exon(s)
var pri_exon = [], n_pri_exons;
{
	var cnt = [], max = 0;
	for (var e = 0; e < gcnt.length; ++e) {
		if (gcnt[e] != null) {
			var c = 0, h = gcnt[e];
			for (var x in h) ++c;
			cnt[e] = c;
			max = max > c? max : c;
		} else cnt[e] = 0;
	}
	var pri_list = [];
	for (var e = 0; e < cnt.length; ++e) {
		if (cnt[e] == max) pri_list.push(e + 1);
		pri_exon[e] = cnt[e] == max? 1 : 0;
	}
	warn("List of primary exons: ["+pri_list.join(",")+"]");
	n_pri_exons = pri_list.length;
}

// convert strings to integers (for performance)
var ghash = {}, glist = [], chash = {}, clist = [], elist = [];
for (var i = 0; i < list.length; ++i) {
	var g = glist.length, c = clist.length;
	if (ghash[list[i][1]] == null) {
		glist.push(list[i][1]);
		ghash[list[i][1]] = g;
	}
	if (chash[list[i][0]] == null) {
		clist.push(list[i][0]);
		chash[list[i][0]] = c;
	}
	if (elist[g] == null) elist[g] = {};
	elist[g][list[i][2]] = true;
}

// change to the per-exon representation
var exons = [];
for (var i = 0; i < list.length; ++i) {
	var li = list[i];
	if (exons[li[2]] == null) exons[li[2]] = [];
	exons[li[2]].push([chash[li[0]], ghash[li[1]], li[5] + li[6]]);
}

// initialize genotype scores
var pair = [];
for (var i = 0; i < glist.length; ++i) {
	pair[i] = [];
	for (var j = 0; j <= i; ++j)
		pair[i][j] = 0;
}

// type each exon
for (var e = 0; e < exons.length; ++e) {
//for (var e = 1; e < 3; ++e) {
	if (exons[e] == null) continue;
	var ee = exons[e];
	// find good contigs and alleles that have good matches
	var ch = {}, gh = {};
	for (var i = 0; i < ee.length; ++i) {
		if (elist[ee[i][1]][e] != null)
			ch[ee[i][0]] = true, gh[ee[i][1]] = true;
	}
	var ca = [], ga = [];
	for (var c in ch) ca.push(parseInt(c));
	for (var g in gh) ga.push(parseInt(g));
	var named_ca = [];
	for (var i = 0; i < ca.length; ++i) named_ca.push(clist[ca[i]]);
	warn("Processing exon "+(e+1)+" (" +ga.length+ " genes; " +ca.length+ " contigs: [" +named_ca.join(",")+ "])...");
	// convert representation again
	var sc = [];
	for (var i = 0; i < ee.length; ++i) {
		var c = ee[i][0], g = ee[i][1];
		if (sc[g] == null) sc[g] = [];
		if (sc[g][c] == null) sc[g][c] = 0xffff;
		sc[g][c] = sc[g][c] < ee[i][2]? sc[g][c] : ee[i][2];
	}
	// set unmapped entries to high mismatch
	for (var i = 0; i < ga.length; ++i)
		for (var j = 0; j < ca.length; ++j) {
			var g = ga[i], c = ca[j];
			if (sc[g][c] == null) sc[g][c] = 0xff;
		}
	// drop mismapped contigs
	var dropped = [];
	for (var c = 0; c < ca.length; ++c) {
		var min = 0x7fffffff, cc = ca[c];
		for (var g = 0; g < ga.length; ++g) {
			var gg = ga[g];
			min = min < sc[gg][cc]? min : sc[gg][cc];
		}
		dropped[cc] = min > drop_thres? true : false;
		if (dropped[cc]) warn("Dropped contig " +clist[cc]+ " due to high divergence to all genes (minNM=" +min+ ")");
	}
	// fill the pair array
	var min_nm = 0xffff;
	for (var i = 0; i < ga.length; ++i) {
		var gi = ga[i], g1 = sc[gi];
		for (var j = i; j < ga.length; ++j) {
			var gj = ga[j], g2 = sc[gj], m = 0;
			for (var k = 0; k < ca.length; ++k) {
				c = ca[k];
				if (!dropped[c])
					m += g1[c] < g2[c]? g1[c] : g2[c];
			}
			var x = m<<20 | 1<<6 | pri_exon[e];
			if (gi < gj) pair[gj][gi] += x;
			else pair[gi][gj] += x;
			min_nm = min_nm < m? min_nm : m;
		}
	}
}

// genotyping
var min_nm = 0x7fffffff;
for (var i = 0; i < glist.length; ++i)
	for (var j = 0; j <= i; ++j)
		if ((pair[i][j]&63) == n_pri_exons)
			min_nm = min_nm < pair[i][j]>>20? min_nm : pair[i][j]>>20;

var out = [];
for (var i = 0; i < glist.length; ++i)
	for (var j = 0; j <= i; ++j)
		if ((pair[i][j]&63) == n_pri_exons && pair[i][j]>>20 <= min_nm + fuzzy)
			out.push([pair[i][j]>>20, pair[i][j]>>6&63, i, j]);

out.sort(function(a, b) { return a[0]!=b[0]? a[0]-b[0] : b[1]-a[1]});
for (var i = 0; i < out.length; ++i)
	print(glist[out[i][2]], glist[out[i][3]], out[i][0], out[i][1]);
