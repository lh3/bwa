/*****************************************************************
 * The K8 Javascript interpreter is required to run this script. *
 *                                                               *
 *   Source code: https://github.com/attractivechaos/k8          *
 *   Binary: http://sourceforge.net/projects/lh3/files/k8/       *
 *****************************************************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/************************
 * Command line parsing *
 ************************/

var ver = "r19";
var c, thres_len = 50, thres_ratio = .8, thres_nm = 5, thres_frac = .33, dbg = false;

// parse command line options
while ((c = getopt(arguments, "vdl:n:f:")) != null) {
	if (c == 'l') thres_len = parseInt(getopt.arg);
	else if (c == 'n') thres_nm = parseInt(getopt.arg);
	else if (c == 'd') dbg = true;
	else if (c == 'f') thres_frac = parseFloat(getopt.arg);
	else if (c == 'v') { print(ver); exit(0); }
}
if (arguments.length == getopt.ind) {
	print("");
	print("Usage:   k8 typeHLA.js [options] <exon-to-contig.sam>\n");
	print("Options: -n INT     drop a contig if the edit distance to the closest gene is >INT ["+thres_nm+"]");
	print("         -l INT     drop a contig if its match too short ["+thres_len+"]");
	print("         -f FLOAT   drop inconsistent contigs if their length <FLOAT fraction of total length ["+thres_ratio.toFixed(2)+"]");
	print("         -d         output extra info for debugging");
	print("         -v         show version number");
	print("");
	print("Note: The output is TAB delimited with each GT line consisting of allele1, allele2,");
	print("      #mismatches/gaps on primary exons, #mismatches/gaps on other exons and #exons");
	print("      used in typing. If unusure, use the first GT line as the final genotype.\n");
	exit(1);
}

/*********************************
 * Read gene-to-contig alignment *
 *********************************/

var file = new File(arguments[getopt.ind]);
var buf = new Bytes();
var re_cigar = /(\d+)([MIDSH])/g;

var len = {}, list = [], gcnt = [];
while (file.readline(buf) >= 0) {
	var m, mm, line = buf.toString();
	var t = line.split("\t");
	var flag = parseInt(t[1]);
	// SAM header
	if (t[0].charAt(0) == '@') {
		if (t[0] == '@SQ' && (m = /LN:(\d+)/.exec(line)) != null && (mm = /SN:(\S+)/.exec(line)) != null)
			len[mm[1]] = parseInt(m[1]);
		continue;
	}
	// parse gene name and exon number
	var gene = null, exon = null;
	if ((m = /^(HLA-[^\s_]+)_(\d+)/.exec(t[0])) != null) {
		gene = m[1], exon = parseInt(m[2]) - 1;
		if (gcnt[exon] == null) gcnt[exon] = {};
		gcnt[exon][gene] = true;
	}
	if (gene == null || exon == null || t[2] == '*') continue;
	// parse clipping and aligned length
	var x = 0, ts = parseInt(t[3]) - 1, te = ts, clip = [0, 0];
	while ((m = re_cigar.exec(t[5])) != null) {
		var l = parseInt(m[1]);
		if (m[2] == 'M') x += l, te += l;
		else if (m[2] == 'I') x += l;
		else if (m[2] == 'D') te += l;
		else if (m[2] == 'S' || m[2] == 'H') clip[x==0?0:1] = l;
	}
	var tl = len[t[2]];
	var left  = ts < clip[0]? ts : clip[0];
	var right = tl - te < clip[1]? tl - te : clip[1];
	var qs, qe, ql = clip[0] + x + clip[1];
	if (flag & 16) qs = clip[1], qe = ql - clip[0];
	else qs = clip[0], qe = ql - clip[1];
	var nm = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : 0;
	list.push([t[2], gene, exon, ts, te, nm, left + right, qs, qe, ql]); // left+right should be 0 given a prefix-suffix alignment
}

buf.destroy();
file.close();

/**************************************
 * Prepare data structures for typing *
 **************************************/

// identify the primary exons, the exons associated with most genes
var pri_exon = [], n_pri_exons;
{
	var cnt = [], max = 0;
	// count the number of genes per exon and track the max
	for (var e = 0; e < gcnt.length; ++e) {
		if (gcnt[e] != null) {
			var c = 0, h = gcnt[e];
			for (var x in h) ++c;
			cnt[e] = c;
			max = max > c? max : c;
		} else cnt[e] = 0;
	}
	warn("- Number of genes for each exon: [" +cnt.join(",") + "]");
	// find primary exons
	var pri_list = [];
	for (var e = 0; e < cnt.length; ++e) {
		if (cnt[e] == max) pri_list.push(e + 1);
		pri_exon[e] = cnt[e] == max? 1 : 0;
	}
	warn("- List of primary exon(s): ["+pri_list.join(",")+"]");
	n_pri_exons = pri_list.length;
}

// convert strings to integers (for performance)
var ghash = {}, glist = [], chash = {}, clist = [], elist = [];
for (var i = 0; i < list.length; ++i) {
	if (ghash[list[i][1]] == null) {
		ghash[list[i][1]] = glist.length;
		glist.push(list[i][1]);
	}
	if (chash[list[i][0]] == null) {
		chash[list[i][0]] = clist.length;
		clist.push(list[i][0]);
	}
	var g = ghash[list[i][1]];
	if (elist[g] == null) elist[g] = {};
	elist[g][list[i][2]] = true;
}

// extract the 3rd and 4th digits
var gsub = [], gsuf = [];
for (var i = 0; i < glist.length; ++i) {
	var m = /^HLA-[^*\s]+\*\d+:(\d+).*([A-Z]?)$/.exec(glist[i]);
	gsub[i] = parseInt(m[1]);
	gsuf[i] = /[A-Z]$/.test(glist[i])? 1 : 0;
}

/*************************************************
 * Collect genes with perfect matches on primary *
 *************************************************/

// collect exons with fully covered by perfect match(es)
var perf_exons = [];

function push_perf_exons(matches, last)
{
	matches.sort(function(a, b) { return a[0]-b[0]; });
	var cov = 0, start = 0, end = 0;
	for (var i = 0; i < matches.length; ++i) {
		if (matches[i][3] > 0) continue;
		if (matches[i][0] <= end)
			end = end > matches[i][1]? end : matches[i][1];
		else cov += end - start, start = matches[i][0], end = matches[i][1];
	}
	cov += end - start;
	if (matches[0][2] == cov) {
		if (perf_exons[last[1]] == null) perf_exons[last[1]] = [];
		//print(last[0], last[1], ghash[last[0]]);
		perf_exons[last[1]].push(ghash[last[0]]);
	}
}

var last = [null, -1], matches = [];
for (var i = 0; i < list.length; ++i) {
	var li = list[i];
	if (last[0] != li[1] || last[1] != li[2]) {
		if (matches.length) push_perf_exons(matches, last);
		matches = [];
		last = [li[1], li[2]];
	}
	matches.push([li[7], li[8], li[9], li[5]+li[6]]);
}
if (matches.length) push_perf_exons(matches, last);

// for each gene, count how many primary exons are perfect
var pg_aux_cnt = {};
for (var e = 0; e < perf_exons.length; ++e) {
	if (!pri_exon[e]) continue;
	var pe = perf_exons[e];
	var n = pe? pe.length : 0;
	for (var i = 0; i < n; ++i) {
		var g = pe[i];
		if (pg_aux_cnt[g] == null) pg_aux_cnt[g] = 1;
		else ++pg_aux_cnt[g];
	}
}

// find genes with perfect matches on the primary exons
var perf_genes = [];
for (var g in pg_aux_cnt)
	if (pg_aux_cnt[g] == n_pri_exons)
		perf_genes.push(parseInt(g));
warn("- Found " +perf_genes.length+ " genes fully covered by perfect matches on the primary exon(s)");

var h_perf_genes = {};
for (var i = 0; i < perf_genes.length; ++i) {
	if (dbg) print("PG", glist[perf_genes[i]]);
	h_perf_genes[perf_genes[i]] = true;
}

/*******************
 * Filter hit list *
 *******************/

// reorganize hits to exons
function list2exons(list, flt_flag, perf_hash)
{
	var exons = [];
	for (var i = 0; i < list.length; ++i) {
		var li = list[i], c = chash[li[0]], g = ghash[li[1]];
		if (flt_flag != null && flt_flag[c] == 1) continue;
		if (perf_hash != null && !perf_hash[g]) continue;
		if (exons[li[2]] == null) exons[li[2]] = [];
		exons[li[2]].push([c, g, li[5] + li[6], li[4] - li[3]]);
	}
	return exons;
}

var exons = list2exons(list), flt_flag = [], ovlp_len = [];
for (var c = 0; c < clist.length; ++c) flt_flag[c] = ovlp_len[c] = 0;
for (var e = 0; e < exons.length; ++e) {
	if (!pri_exon[e]) continue;
	var ee = exons[e];
	var max_len = [];
	for (var c = 0; c < clist.length; ++c) max_len[c] = 0;
	for (var i = 0; i < ee.length; ++i) {
		var l = ee[i][3] - ee[i][2];
		if (l < 1) l = 1;
		if (max_len[ee[i][0]] < l) max_len[ee[i][0]] = l;
	}
	for (var c = 0; c < clist.length; ++c) ovlp_len[c] += max_len[c];
	for (var i = 0; i < ee.length; ++i)
		flt_flag[ee[i][0]] |= (!h_perf_genes[ee[i][1]] || ee[i][2])? 1 : 1<<1;	
}

var l_cons = 0, l_incons = 0;
for (var c = 0; c < clist.length; ++c)
	if (flt_flag[c]&2) l_cons += ovlp_len[c];
	else if (flt_flag[c] == 1) l_incons += ovlp_len[c];

warn("- Total length of contigs consistent/inconsistent with perfect genes: " +l_cons+ "/" +l_incons);
var attempt_perf = (l_incons/(l_cons+l_incons) < thres_frac);

/********************************
 * Core function for genotyping *
 ********************************/

function type_gene(perf_mode)
{
	if (perf_mode) {
		var flt_list = [];
		for (var c = 0; c < clist.length; ++c)
			if (flt_flag[c] == 1) flt_list.push(clist[c]);
		warn("  - Filtered " +flt_list.length+ " inconsistent contig(s): [" +flt_list.join(",")+ "]");
		exons = list2exons(list, flt_flag, h_perf_genes);
	} else exons = list2exons(list);

	/***********************
	 * Score each genotype *
	 ***********************/

	// initialize genotype scores
	var pair = [];
	for (var i = 0; i < glist.length; ++i) {
		pair[i] = [];
		for (var j = 0; j <= i; ++j)
			pair[i][j] = 0;
	}

	// these two arrays are used to output debugging information
	var score = [], ctg = [];

	function type_exon(e, gt_list)
	{
		function update_pair(x, m, is_pri)
		{
			var y, z;
			y = (x>>14&0xff) + m < 0xff? (x>>14&0xff) + m : 0xff;
			if (is_pri) z = (x>>22) + m < 0xff? (x>>22) + m : 0xff;
			else z = x>>22;
			return z<<22 | y<<14 | ((x&0x3fff) + (1<<6|is_pri));
		}

		score[e] = []; ctg[e] = [];
		if (exons[e] == null) return;
		var ee = exons[e], is_pri = pri_exon[e]? 1 : 0;
		// find contigs and genes associated with the current exon
		var ch = {}, gh = {};
		for (var i = 0; i < ee.length; ++i)
			if (elist[ee[i][1]][e] != null)
				ch[ee[i][0]] = true, gh[ee[i][1]] = true;
		var ga = [], ca = ctg[e];
		for (var c in ch) ca.push(parseInt(c));
		for (var g in gh) ga.push(parseInt(g));
		var named_ca = [];
		for (var i = 0; i < ca.length; ++i) named_ca.push(clist[ca[i]]);
		warn("    - Processing exon "+(e+1)+" (" +ga.length+ " genes; " +ca.length+ " contigs: [" +named_ca.join(", ")+ "])...");
		// set unmapped entries to high mismatch
		var sc = score[e];
		for (var k = 0; k < ga.length; ++k) {
			var g = ga[k];
			if (sc[g] == null) sc[g] = [];
			for (var i = 0; i < ca.length; ++i)
				sc[g][ca[i]] = 0xff;
		}
		// convert representation again and compute max_len[]
		var max_len = [];
		for (var i = 0; i < ee.length; ++i) {
			var c = ee[i][0], g = ee[i][1];
			if (gh[g] == null || ch[c] == null) continue;
			sc[g][c] = sc[g][c] < ee[i][2]? sc[g][c] : ee[i][2];
			if (max_len[c] == null) max_len[c] = 0;
			max_len[c] = max_len[c] > ee[i][3]? max_len[c] : ee[i][3];
		}
		// drop mismapped contigs
		var max_max_len = 0;
		for (var k = 0; k < ca.length; ++k)
			max_max_len = max_max_len > max_len[ca[k]]? max_max_len : max_len[ca[k]];
		var dropped = [];
		for (var k = 0; k < ca.length; ++k) {
			var min = 0x7fffffff, c = ca[k];
			for (var i = 0; i < ga.length; ++i) {
				var g = ga[i];
				min = min < sc[g][c]? min : sc[g][c];
			}
			dropped[c] = min > thres_nm? true : false;
			if (max_len[c] < thres_len && max_len[c] < thres_ratio * max_max_len) dropped[c] = true;
			if (dropped[c]) warn("      . Dropped low-quality contig " +clist[c]+ " (minNM=" +min+ "; maxLen=" +max_len[c]+ ")");
		}
		// fill the pair array
		if (gt_list == null) {
			for (var i = 0; i < ga.length; ++i) {
				var m = 0, gi = ga[i], g1 = sc[gi];
				// homozygous
				for (var k = 0; k < ca.length; ++k) {
					var c = ca[k];
					if (!dropped[c]) m += g1[c];
				}
				pair[gi][gi] = update_pair(pair[gi][gi], m, is_pri);
				// heterozygous
				for (var j = i + 1; j < ga.length; ++j) {
					var gj = ga[j], g2 = sc[gj], m = 0, a = [0, 0];
					for (var k = 0; k < ca.length; ++k) {
						var c = ca[k];
						if (!dropped[c]) {
							m += g1[c] < g2[c]? g1[c] : g2[c];
							++a[g1[c]<g2[c]? 0:1];
						}
					}
					if (a[0] == 0 || a[1] == 0) m = 0xff; // if all contigs are assigned to one gene, it is not good
					if (gi < gj) pair[gj][gi] = update_pair(pair[gj][gi], m, is_pri);
					else pair[gi][gj] = update_pair(pair[gi][gj], m, is_pri);
				}
			}
		} else {
			var tmp_pairs = [], min = 0xff;
			for (var i = 0; i < gt_list.length; ++i) {
				var gt = gt_list[i], m = 0;
				var g1 = sc[gt[0]], g2 = sc[gt[1]], a = [0, 0];
				if (g1 == null || g2 == null) continue;
				if (gt[0] == gt[1]) {
					for (var k = 0; k < ca.length; ++k) {
						var c = ca[k];
						if (!dropped[c]) m += g1[c];
					}
				} else {
					var a = [0, 0];
					for (k = 0; k < ca.length; ++k) {
						var c = ca[k];
						if (!dropped[c]) {
							m += g1[c] < g2[c]? g1[c] : g2[c];
							++a[g1[c]<g2[c]? 0:1];
						}
					}
					if (a[0] == 0 || a[1] == 0) m = 0xff;
				}
				tmp_pairs.push([gt[0], gt[1], m]);
				min = min < m? min : m;
			}
			if (min < 0xff) {
				for (var i = 0; i < tmp_pairs.length; ++i) {
					var t = tmp_pairs[i];
					pair[t[0]][t[1]] = update_pair(pair[t[0]][t[1]], t[2], is_pri);
				}
			} else warn("      . Skipped exon " +(e+1)+ " as the assembly may be incomplete");
		}
	}

	// type primary exons
	warn("  - Processing primary exon(s)...");
	for (var e = 0; e < exons.length; ++e)
		if (pri_exon[e]) type_exon(e);

	// generate the list of best genotypes on primary exons
	var min_nm_pri = 0x7fffffff;
	for (var i = 0; i < glist.length; ++i)
		for (var j = 0; j <= i; ++j)
			if ((pair[i][j]&63) == n_pri_exons)
				min_nm_pri = min_nm_pri < pair[i][j]>>22? min_nm_pri : pair[i][j]>>22;

	var gt_list = [];
	for (var i = 0; i < glist.length; ++i)
		for (var j = 0; j <= i; ++j)
			if ((pair[i][j]&63) == n_pri_exons && pair[i][j]>>22 == min_nm_pri)
				gt_list.push([i, j]);

	warn("  - Collected " +gt_list.length+ " top genotypes on the primary exon(s); minimal edit distance: " +min_nm_pri);

	// type other exons
	warn("  - Processing other exon(s)...");
	for (var e = 0; e < exons.length; ++e)
		if (!pri_exon[e]) type_exon(e, gt_list);

	/*****************************
	 * Choose the best genotypes *
	 *****************************/

	// genotyping
	var min_nm = 0x7fffffff;
	for (var i = 0; i < glist.length; ++i)
		for (var j = 0; j <= i; ++j)
			if ((pair[i][j]&63) == n_pri_exons)
				min_nm = min_nm < pair[i][j]>>14? min_nm : pair[i][j]>>14;

	var out = [];
	for (var i = 0; i < glist.length; ++i)
		for (var j = 0; j <= i; ++j)
			if ((pair[i][j]&63) == n_pri_exons && pair[i][j]>>14 <= min_nm + 1)
				out.push([pair[i][j]>>14, pair[i][j]>>6&0xff, i, j, (gsuf[i] + gsuf[j])<<16|(gsub[i] + gsub[j])]);

	out.sort(function(a, b) { return a[0]!=b[0]? a[0]-b[0] : a[1]!=b[1]? b[1]-a[1] : a[4]!=b[4]? a[4]-b[4] : a[2]!=b[2]? a[2]-b[2] : a[3]-b[3]});

	return out;
}

/**********************
 * Perform genotyping *
 **********************/

warn("- Typing in the imperfect mode...");
var rst = type_gene(false);
if (attempt_perf) {
	warn("- Typing in the perfect mode...");
	var rst_perf = type_gene(true);
	warn("- Imperfect vs perfect mode: [" +(rst[0][0]>>8&0xff)+ "," +(rst[0][0]&0xff)+ "] vs [" +(rst_perf[0][0]>>8&0xff)+ "," +(rst_perf[0][0]&0xff)+ "]");
	if (rst_perf[0][0] < rst[0][0]) {
		warn("- Chose the result from the perfect mode");
		rst = rst_perf;
	} else warn("- Chose the result from the imperfect mode");
} else warn("- Perfect mode is not attempted");

/**********
 * Output *
 **********/

for (var i = 0; i < rst.length; ++i)
	print("GT", glist[rst[i][3]], glist[rst[i][2]], rst[i][0]>>8&0xff, rst[i][0]&0xff, rst[i][1]);
