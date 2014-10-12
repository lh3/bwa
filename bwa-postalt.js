/*****************************************************************
 * The K8 Javascript interpreter is required to run this script. *
 *                                                               *
 *   Source code: https://github.com/attractivechaos/k8          *
 *   Binary: http://sourceforge.net/projects/lh3/files/k8/       *
 *                                                               *
 * Data file used for generating GRCh38 ALT alignments:          *
 *                                                               *
 *   http://sourceforge.net/projects/bio-bwa/files/              *
 *****************************************************************/

/******************
 *** From k8.js ***
 ******************/

// Parse command-line options. A BSD getopt() clone in javascript.
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

// print an object in a format similar to JSON. For debugging only.
function obj2str(o)
{
	if (typeof(o) != 'object') {
		return o.toString();
	} else if (o == null) {
		return "null";
	} else if (Array.isArray(o)) {
		var s = "[";
		for (var i = 0; i < o.length; ++i) {
			if (i) s += ',';
			s += obj2str(o[i]);
		}
		return s + "]";
	} else {
		var i = 0, s = "{";
		for (var key in o) {
			if (i++) s += ',';
			s += key + ":";
			s += obj2str(o[key]);
		}
		return s + "}";
	}
}

// reverse a string
Bytes.prototype.reverse = function()
{
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[i];
		this[i] = this[this.length - i - 1];
		this[this.length - i - 1] = tmp;
	}
}

// reverse complement a DNA string
Bytes.prototype.revcomp = function()
{
	if (Bytes.rctab == null) {
		var s1 = 'WSATUGCYRKMBDHVNwsatugcyrkmbdhvn';
		var s2 = 'WSTAACGRYMKVHDBNwstaacgrymkvhdbn';
		Bytes.rctab = [];
		for (var i = 0; i < 256; ++i) Bytes.rctab[i] = 0;
		for (var i = 0; i < s1.length; ++i)
			Bytes.rctab[s1.charCodeAt(i)] = s2.charCodeAt(i);
	}
	for (var i = 0; i < this.length>>1; ++i) {
		var tmp = this[this.length - i - 1];
		this[this.length - i - 1] = Bytes.rctab[this[i]];
		this[i] = Bytes.rctab[tmp];
	}
	if (this.length>>1)
		this[this.length>>1] = Bytes.rctab[this[this.length>>1]];
}

// create index for a list of intervals for fast interval queries; ported from bedidx.c in samtools
function intv_ovlp(intv, bits)
{
	if (typeof bits == "undefined") bits = 13;
	intv.sort(function(a,b) {return a[0]-b[0];});
	// create the index
	var idx = [], max = 0;
	for (var i = 0; i < intv.length; ++i) {
		var b = intv[i][0]>>bits;
		var e = (intv[i][1]-1)>>bits;
		if (b != e) {
			for (var j = b; j <= e; ++j)
				if (idx[j] == null) idx[j] = i;
		} else if (idx[b] == null) idx[b] = i;
		max = max > e? max : e;
	}
	// closure
	return function(_b, _e) {
		var x = _b >> bits;
		if (x > max) return [];
		var off = idx[x];
		if (off == null) {
			var i;
			for (i = ((_e - 1) >> bits) - 1; i >= 0; --i)
				if (idx[i] != null) break;
			off = i < 0? 0 : idx[i];
		}
		var ovlp = [];
		for (var i = off; i < intv.length && intv[i][0] < _e; ++i)
			if (intv[i][1] > _b) ovlp.push(intv[i]);
		return ovlp;
	}
}

var re_cigar = /(\d+)([MIDSHN])/g;

/******************************
 *** Generate ALT alignment ***
 ******************************/

// given a pos on ALT and the ALT-to-REF CIGAR, find the pos on REF
function cigar2pos(cigar, pos)
{
	var x = 0, y = 0;
	for (var i = 0; i < cigar.length; ++i) {
		var op = cigar[i][0], len = cigar[i][1];
		if (op == 'M') {
			if (y <= pos && pos < y + len)
				return x + (pos - y);
			x += len, y += len;
		} else if (op == 'D') {
			x += len;
		} else if (op == 'I') {
			if (y <= pos && pos < y + len)
				return x;
			y += len;
		} else if (op == 'S' || op == 'H') {
			if (y <= pos && pos < y + len)
				return -1;
			y += len;
		}
	}
	return -1;
}

// Parse a hit. $s is an array that looks something like ["chr1", "+12345", "100M", 5]
// Return an object keeping various information about the alignment.
function parse_hit(s, opt)
{
	var h = {};
	h.ctg = s[0];
	h.start = parseInt(s[1].substr(1)) - 1;
	h.rev = (s[1].charAt(0) == '-');
	h.cigar = s[2];
	h.NM = parseInt(s[3]);
	h.hard = false;
	var m, l_ins, n_ins, l_del, n_del, l_match, l_skip, l_clip;
	l_ins = l_del = n_ins = n_del = l_match = l_skip = l_clip = 0;
	while ((m = re_cigar.exec(h.cigar)) != null) {
		var l = parseInt(m[1]);
		if (m[2] == 'M') l_match += l;
		else if (m[2] == 'D') ++n_del, l_del += l;
		else if (m[2] == 'I') ++n_ins, l_ins += l;
		else if (m[2] == 'N') l_skip += l;
		else if (m[2] == 'H' || m[2] == 'S') {
			l_clip += l;
			if (m[2] == 'H') h.hard = true;
		}
	}
	h.end = h.start + l_match + l_del + l_skip;
	h.NM = h.NM > l_del + l_ins? h.NM : l_del + l_ins;
	h.score = Math.floor((opt.a * l_match - (opt.a + opt.b) * (h.NM - l_del - l_ins) - opt.o * (n_del + n_ins) - opt.e * (l_del + l_ins)) / opt.a + .499);
	h.l_query = l_match + l_ins + l_clip;
	return h;
}

function bwa_postalt(args)
{
	var c, opt = { a:1, b:4, o:6, e:1, verbose:3, show_pri:false, recover_mapq:true, min_mapq:10, min_sc:90, max_nm_sc:100, show_ev:false, wei1:false };

	while ((c = getopt(args, '1pqev:')) != null) {
		if (c == 'v') opt.verbose = parseInt(getopt.arg);
		else if (c == 'p') opt.show_pri = true;
		else if (c == 'q') opt.recover_maq = false;
		else if (c == 'e') opt.show_ev = true;
		else if (c == '1') opt.wei1 = true;
	}

	if (args.length == getopt.ind) {
		print("");
		print("Usage:   k8 bwa-postalt.js [-p] <alt.sam> [aln.sam]\n");
		print("Options: -p    output lifted non-ALT hit in a SAM line (for ALT-unware alignments)");
		print("         -q    don't recover mapQ for non-ALTs hit overlapping lifted ALT");
		print("");
		print("Note: This script inspects the XA tag, lifts the mapping positions of ALT hits to");
		print("      the primary assembly, groups them and then estimates mapQ across groups. If");
		print("      a non-ALT hit overlaps a lifted ALT hit, its mapping quality is set to the");
		print("      smaller between its original mapQ and the adjusted mapQ of the ALT hit. If");
		print("      multiple ALT hits are lifted to the same position, they will yield new SAM");
		print("      lines with the same mapQ.");
		print("");
		exit(1);
	}

	var aux = new Bytes(); // used for reverse and reverse complement
	var buf = new Bytes();

	// read ALT-to-REF alignment
	var intv_alt = {}, intv_pri = {}, idx_un = {};
	var file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line.charAt(0) == '@') continue;
		var t = line.split("\t");
		if (t.length < 11) continue; // incomplete lines
		var pos = parseInt(t[3]) - 1;
		var flag = parseInt(t[1]);
		if ((flag&4) || t[2] == '*') {
			idx_un[t[0]] = true;
			continue;
		}
		var m, cigar = [], l_qaln = 0, l_tlen = 0, l_qclip = 0;
		while ((m = re_cigar.exec(t[5])) != null) {
			var l = parseInt(m[1]);
			cigar.push([m[2] != 'H'? m[2] : 'S', l]); // convert hard clip to soft clip
			if (m[2] == 'M') l_qaln += l, l_tlen += l;
			else if (m[2] == 'I') l_qaln += l;
			else if (m[2] == 'S' || m[2] == 'H') l_qclip += l;
			else if (m[2] == 'D' || m[2] == 'N') l_tlen += l;
		}
		var j = flag&16? cigar.length-1 : 0;
		var start = cigar[j][0] == 'S'? cigar[j][1] : 0;
		if (intv_alt[t[0]] == null) intv_alt[t[0]] = [];
		intv_alt[t[0]].push([start, start + l_qaln, l_qaln + l_qclip, t[2], flag&16? true : false, pos - 1, cigar, pos + l_tlen]);
		if (intv_pri[t[2]] == null) intv_pri[t[2]] = [];
		intv_pri[t[2]].push([pos, pos + l_tlen, t[0]]);
	}
	file.close();
	var idx_alt = {}, idx_pri = {};
	for (var ctg in intv_alt)
		idx_alt[ctg] = intv_ovlp(intv_alt[ctg]);
	for (var ctg in intv_pri)
		idx_pri[ctg] = intv_ovlp(intv_pri[ctg]);

	// initialize the list of ALT contigs
	var weight_alt = [];
	for (var ctg in idx_alt)
		weight_alt[ctg] = [0, 0, 0, intv_alt[ctg][0][3], intv_alt[ctg][0][5], intv_alt[ctg][0][7]];
	for (var ctg in idx_un)
		weight_alt[ctg] = [0, 0, 0, '~', 0, 0];

	// process SAM
	var buf2 = [];
	file = args.length - getopt.ind >= 2? new File(args[getopt.ind+1]) : new File();
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if (line.charAt(0) == '@') { // print and then skip the header line
			print(line);
			continue;
		}

		var t = line.split("\t");
		t[1] = parseInt(t[1]); t[3] = parseInt(t[3]); t[4] = parseInt(t[4]);

		// print bufferred reads
		if (buf2.length && (buf2[0][0] != t[0] || (buf2[0][1]&0xc0) != (t[1]&0xc0))) {
			for (var i = 0; i < buf2.length; ++i)
				print(buf2[i].join("\t"));
			buf2 = [];
		}

		// skip unmapped lines
		if (t[1]&4) {
			buf2.push(t);
			continue;
		}

		// parse the reported hit
		var NM = (m = /\tNM:i:(\d+)/.exec(line)) == null? '0' : m[1];
		var flag = t[1];
		var h = parse_hit([t[2], ((flag&16)?'-':'+') + t[3], t[5], NM], opt);
		if (h.hard) { // the following does not work with hard clipped alignments
			buf2.push(t);
			continue;
		}
		var hits = [h];

		// parse hits in the XA tag
		if ((m = /\tXA:Z:(\S+)/.exec(line)) != null) {
			var XA_strs = m[1].split(";");
			for (var i = 0; i < XA_strs.length; ++i)
				if (XA_strs[i] != '') // as the last symbol in an XA tag is ";", the last split is an empty string
					hits.push(parse_hit(XA_strs[i].split(","), opt));
		}

		// check if there are ALT hits
		var has_alt = false;
		for (var i = 0; i < hits.length; ++i)
			if (weight_alt[hits[i].ctg] != null) {
				has_alt = true;
				break;
			}
		if (!has_alt) {
			buf2.push(t);
			continue;
		}

		// lift mapping positions to the primary assembly
		var n_rpt_lifted = 0, rpt_lifted = null;
		for (var i = 0; i < hits.length; ++i) {
			var a, h = hits[i];

			if (idx_alt[h.ctg] == null || (a = idx_alt[h.ctg](h.start, h.end)) == null || a.length == 0)
				continue;

			// find the approximate position on the primary assembly
			var lifted = [];
			for (var j = 0; j < a.length; ++j) {
				var s, e;
				if (!a[j][4]) { // ALT is mapped to the forward strand of the primary assembly
					s = cigar2pos(a[j][6], h.start);
					e = cigar2pos(a[j][6], h.end - 1) + 1;
				} else {
					s = cigar2pos(a[j][6], a[j][2] - h.end);
					e = cigar2pos(a[j][6], a[j][2] - h.start - 1) + 1;
				}
				if (s < 0 || e < 0) continue; // read is mapped to clippings in the ALT-to-chr alignment
				s += a[j][5]; e += a[j][5];
				lifted.push([a[j][3], (h.rev!=a[j][4]), s, e]);
				if (i == 0) ++n_rpt_lifted;
			}
			if (i == 0 && n_rpt_lifted == 1) rpt_lifted = lifted[0].slice(0);
			if (lifted.length) hits[i].lifted = lifted;
		}

		// prepare for hits grouping
		for (var i = 0; i < hits.length; ++i) { // set keys for sorting
			if (hits[i].lifted != null) // TODO: only the first element in lifted[] is used
				hits[i].pctg = hits[i].lifted[0][0], hits[i].pstart = hits[i].lifted[0][2], hits[i].pend = hits[i].lifted[0][3];
			else hits[i].pctg = hits[i].ctg, hits[i].pstart = hits[i].start, hits[i].pend = hits[i].end;
			hits[i].i = i; // keep the original index
		}

		// group hits based on the lifted positions on non-ALT sequences
		if (hits.length > 1) {
			hits.sort(function(a,b) { return a.pctg != b.pctg? (a.pctg < b.pctg? -1 : 1) : a.pstart - b.pstart });
			var last_chr = null, end = 0, g = -1;
			for (var i = 0; i < hits.length; ++i) {
				if (last_chr != hits[i].pctg) ++g, last_chr = hits[i].pctg, end = 0;
				else if (hits[i].pstart >= end) ++g;
				hits[i].g = g;
				end = end > hits[i].pend? end : hits[i].pend;
			}
		} else hits[0].g = 0;

		// find the index and group id of the reported hit; find the size of the reported group
		var reported_g = null, reported_i = null, n_group0 = 0;
		if (hits.length > 1) {
			for (var i = 0; i < hits.length; ++i)
				if (hits[i].i == 0)
					reported_g = hits[i].g, reported_i = i;
			for (var i = 0; i < hits.length; ++i)
				if (hits[i].g == reported_g)
					++n_group0;
		} else {
			if (weight_alt[hits[0].ctg] == null) { // no need to go through the following if the single hit is non-ALT
				buf2.push(t);
				continue;
			}
			reported_g = reported_i = 0, n_group0 = 1;
		}

		// re-estimate mapping quality if necessary
		var mapQ, ori_mapQ = t[4];
		if (n_group0 > 1) {
			var group_max = [];
			for (var i = 0; i < hits.length; ++i) {
				var g = hits[i].g;
				if (group_max[g] == null || group_max[g][0] < hits[i].score)
					group_max[g] = [hits[i].score, g];
			}
			if (group_max.length > 1)
				group_max.sort(function(x,y) {return y[0]-x[0]});
			if (group_max[0][1] == reported_g) { // the best hit is the hit reported in SAM
				mapQ = group_max.length == 1? 60 : 6 * (group_max[0][0] - group_max[1][0]);
			} else mapQ = 0;
			mapQ = mapQ < 60? mapQ : 60;
			if (idx_alt[t[2]] == null) mapQ = mapQ < ori_mapQ? mapQ : ori_mapQ;
			else mapQ = mapQ > ori_mapQ? mapQ : ori_mapQ;
		} else mapQ = t[4];

		// ALT genotyping
		if (mapQ >= opt.min_mapq && hits[reported_i].score >= opt.min_sc) {
			// collect all overlapping ALT contigs
			var hits2 = [];
			for (var i = 0; i < hits.length; ++i) {
				var h = hits[i];
				if (h.g == reported_g)
					hits2.push([h.pctg, h.pstart, h.pend, h.ctg, h.score, h.NM]);
			}
			var start = hits2[0][1], end = hits2[0][2];
			for (var i = 1; i < hits2.length; ++i)
				end = end > hits2[i][2]? end : hits2[i][2];
			var alts = {};
			for (var i = 0; i < hits2.length; ++i)
				if (weight_alt[hits2[i][3]] != null)
					alts[hits2[i][3]] = [hits2[i][4], hits2[i][5]];
			if (idx_pri[hits2[0][0]] != null) {
				var ovlp = idx_pri[hits2[0][0]](start, end);
				for (var i = 0; i < ovlp.length; ++i) // TODO: requiring reasonable overlap
					if (start <= ovlp[i][0] && ovlp[i][1] <= end && alts[ovlp[i][2]] == null)
						alts[ovlp[i][2]] = [0, 0];
			}

			// add weight to each ALT contig
			var alt_arr = [], max_sc = -1, max_i = -1, sum = 0, min_sc = 1<<30, max_nm = -1;
			for (var ctg in alts)
				alt_arr.push([ctg, alts[ctg][0], 0, alts[ctg][1]]);
			for (var i = 0; i < alt_arr.length; ++i) {
				if (max_sc < alt_arr[i][1])
					max_sc = alt_arr[i][1], max_i = i;
				min_sc = min_sc < alt_arr[i][1]? min_sc : alt_arr[i][1];
				var nm = alt_arr[i][1] > 0? alt_arr[i][3] : opt.max_nm_sc;
				max_nm = max_nm > nm? max_nm : nm;
			}
			if (max_nm > opt.max_nm_sc) max_nm = opt.max_nm_sc;
			if (max_sc > 0 && (alt_arr.length == 1 || min_sc < max_sc)) {
				for (var i = 0; i < alt_arr.length; ++i)
					sum += (alt_arr[i][2] = Math.pow(10, .6 * (alt_arr[i][1] - max_sc)));
				for (var i = 0; i < alt_arr.length; ++i) alt_arr[i][2] /= sum;
				for (var i = 0; i < alt_arr.length; ++i) {
					if (opt.wei1) max_nm = 1;
					var e = [alt_arr[i][0], max_nm, max_nm * alt_arr[max_i][2], max_nm * alt_arr[i][2]];
					var w = weight_alt[e[0]];
					w[0] += e[1], w[1] += e[2], w[2] += e[3];
					if (opt.show_ev) warn(t[0] + '/' + (t[1]>>6&3), e.join("\t"));
				}
			}
		}

		// check if the reported hit overlaps a hit to the primary assembly; if so, don't reduce mapping quality
		if (opt.recover_mapq && n_rpt_lifted == 1 && mapQ > 0) {
			var l = rpt_lifted;
			for (var i = 0; i < buf2.length; ++i) {
				var s = buf2[i];
				if (l[0] != s[2]) continue; // different chr
				if (((s[1]&16) != 0) != l[1]) continue; // different strand
				var start = s[3] - 1, end = start;
				while ((m = re_cigar.exec(t[5])) != null)
					if (m[2] == 'M' || m[2] == 'D' || m[2] == 'N')
						end += parseInt(m[1]);
				if (start < l[3] && l[2] < end) {
					var om = -1;
					for (var j = 11; j < s.length; ++j)
						if ((m = /^om:i:(\d+)/.exec(s[j])) != null)
							om = parseInt(m[1]);
					if (om > 0) s[4] = om;
					s[4] = s[4] < mapQ? s[4] : mapQ;
				}
			}
		}

		// generate lifted_str
		for (var i = 0; i < hits.length; ++i) {
			if (hits[i].lifted && hits[i].lifted.length) {
				var u = '', lifted = hits[i].lifted;
				for (var j = 0; j < lifted.length; ++j)
					u += lifted[j][0] + "," + lifted[j][2] + "," + lifted[j][3] + "," + (lifted[j][1]?'-':'+') + ";";
				hits[i].lifted_str = u;
			}
		}

		// generate reversed quality and reverse-complemented sequence if necessary
		var rs = null, rq = null; // reversed quality and reverse complement sequence
		var need_rev = false;
		for (var i = 0; i < hits.length; ++i) {
			if (hits[i].g != reported_g || i == reported_i) continue;
			if (hits[i].rev != hits[reported_i].rev)
				need_rev = true;
		}
		if (need_rev) { // reverse and reverse complement
			aux.length = 0;
			aux.set(t[9], 0); aux.revcomp(); rs = aux.toString();
			aux.set(t[10],0); aux.reverse(); rq = aux.toString();
		}

		// stage the reported hit
		t[4] = mapQ;
		if (n_group0 > 1) t.push("om:i:"+ori_mapQ);
		if (hits[reported_i].lifted_str) t.push("lt:Z:" + hits[reported_i].lifted_str);
		buf2.push(t);

		// stage the hits generated from the XA tag
		var cnt = 0;
		for (var i = 0; i < hits.length; ++i) {
			if (opt.verbose >= 5) print(obj2str(hits[i]));
			if (hits[i].g != reported_g || i == reported_i) continue;
			if (!opt.show_pri && idx_alt[hits[i].ctg] == null) continue;
			var s = [t[0], 0, hits[i].ctg, hits[i].start+1, mapQ, hits[i].cigar, '*', 0, 0];
			// print sequence/quality and set the rev flag
			if (hits[i].rev == hits[reported_i].rev) {
				s.push(t[9], t[10]);
				s[1] = (flag & 0x10) | 0x800;
			} else {
				s.push(rs, rq);
				s[1] = ((flag & 0x10) ^ 0x10) | 0x800;
			}
			s.push("NM:i:" + hits[i].NM);
			if (hits[i].lifted_str) s.push("lt:Z:" + hits[i].lifted_str);
			buf2.push(s);
		}
	}
	for (var i = 0; i < buf2.length; ++i)
		print(buf2[i].join("\t"));
	file.close();

	buf.destroy();
	aux.destroy();

	// print weight of each contig
	var weight_arr = [];
	for (var ctg in weight_alt) {
		var w = weight_alt[ctg];
		w[0] = w[0].toFixed(4), w[1] = w[1].toFixed(4), w[2] = w[2].toFixed(4);
		weight_arr.push([ctg, w[0], w[1], w[2], w[1] > 0? (w[2]/w[1]).toFixed(3) : '0.000', w[3], w[4], w[5]]);
	}
	weight_arr.sort(function(a,b) {
		return a[5] < b[5]? -1 : a[5] > b[5]? 1 : a[6] != b[6]? a[6] - b[6] : a[0] < b[0]? -1 : a[0] > b[0]? 1 : 0;
	});
	for (var i = 0; i < weight_arr.length; ++i) {
		if (weight_arr[i][5] == '~') weight_arr[i][5] = '*';
		warn(weight_arr[i].join("\t"));
	}
}

bwa_postalt(arguments);
