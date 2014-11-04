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

/*********************
 *** Main function ***
 *********************/

var c, min_qal = 100, drop_thres = 5;

// parse command line options
while ((c = getopt(arguments, "m:d:")) != null) {
	if (c == 'm') min_qal = parseInt(getopt.arg);
	else if (c == 'd') drop_thres = parseInt(getopt.arg);
}
if (arguments.length == getopt.ind) {
	print("");
	print("Usage:   k8 bwa-typeHLA.js [options] <HLA-to-contig.sam>\n");
	print("Options: -d INT     drop a contig if the edit distance to the closest gene is >INT ["+drop_thres+"]");
	print("");
	exit(1);
}

// read alignments
var file = new File(arguments[getopt.ind]);
var buf = new Bytes();
var re_cigar = /(\d+)([MIDSH])/g;

var len = {}, list = [], gcnt = [];
while (file.readline(buf) >= 0) {
	var m, mm, line = buf.toString();
	var t = line.split("\t");
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
	var nm = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : 0;
	list.push([t[2], gene, exon, ts, te, nm, left + right]); // left+right should be 0 given a prefix-suffix alignment
}

buf.destroy();
file.close();

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
	// find primary exons
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

// reorganize hits to exons
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

function update_pair(x, m, is_pri)
{
	var y, z;
	y = (x>>14&0xff) + m < 0xff? (x>>14&0xff) + m : 0xff;
	if (is_pri) z = (x>>22) + m < 0xff? (x>>22) + m : 0xff;
	else z = x>>22;
	return z<<22 | y<<14 | ((x&0x3fff) + (1<<6|is_pri));
}

// type each exon
for (var e = 0; e < exons.length; ++e) {
	if (exons[e] == null) continue;
	var ee = exons[e], is_pri = pri_exon[e]? 1 : 0;
	// find contigs and genes associated with the current exon
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
	warn("Processing exon "+(e+1)+" (" +ga.length+ " genes; " +ca.length+ " contigs: [" +named_ca.join(", ")+ "])...");
	// set unmapped entries to high mismatch
	var sc = [];
	for (var k = 0; k < ga.length; ++k) {
		var g = ga[k];
		if (sc[g] == null) sc[g] = [];
		for (var i = 0; i < ca.length; ++i)
			sc[g][ca[i]] = 0xff;
	}
	// convert representation again
	for (var i = 0; i < ee.length; ++i) {
		var c = ee[i][0], g = ee[i][1];
		sc[g][c] = sc[g][c] < ee[i][2]? sc[g][c] : ee[i][2];
	}
	// drop mismapped contigs
	var dropped = [];
	for (var k = 0; k < ca.length; ++k) {
		var min = 0x7fffffff, c = ca[k];
		for (var i = 0; i < ga.length; ++i) {
			var g = ga[i];
			min = min < sc[g][c]? min : sc[g][c];
		}
		dropped[c] = min > drop_thres? true : false;
		if (dropped[c]) warn("Dropped contig " +clist[c]+ " due to high divergence to all genes (minNM=" +min+ ")");
	}
	// fill the pair array
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
}

// extract the 3rd and 4th digits
var gsub = [], gsuf = [];
for (var i = 0; i < glist.length; ++i) {
	var m = /^HLA-[^*\s]+\*\d+:(\d+).*([A-Z]?)$/.exec(glist[i]);
	gsub[i] = parseInt(m[1]);
	gsuf[i] = /[A-Z]$/.test(glist[i])? 1 : 0;
}

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
for (var i = 0; i < out.length; ++i)
	print(glist[out[i][2]], glist[out[i][3]], out[i][0]>>8&0xff, out[i][0]&0xff, out[i][1]);
