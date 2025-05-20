# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250218
Description:
function position calculator.
"""

class gene_feature_distance_calculator:
    """
    calculate the relative and absolute distances from UTR and CDS endpoints
    args:
        bed_file: path to the BED file
        regions_file: path to the UTR/CDS regions file
        output: path to the output file
    output:
        a bed file with the following columns:
            chr: chromosome
            coord: coordinate
            gene_name: gene name
            refseqID: refseq ID
            rel_location: relative distance from 5'UTR, CDS, or 3'UTR
            utr5_st: start position of 5'UTR
            utr5_end: end position of 5'UTR
            cds_st: start position of CDS
            cds_end: end position of CDS
            utr3_st: start position of 3'UTR
            utr3_end: end position of 3'UTR
            utr5_size: size of 5'UTR
            cds_size: size of CDS
            utr3_size: size of 3'UTR
    """
    def __init__(self, bed_file, regions_file, output):
        """初始化，加载 BED 文件和 UTR/CDS 文件"""
        self.bed_file = bed_file
        self.regions_file = regions_file
        self.refseq_endpts = self.load_refseq_endpoints()
        self.output = output

    def load_refseq_endpoints(self):
        """加载 5'UTR, CDS, 3'UTR 的坐标数据"""
        refseq_endpts = {}
        with open(self.regions_file, 'r') as f:
            for line in f:
                if "\tNA" in line:  # Excludes non-coding RNA (they lack 5'UTR, CDS, and 3'UTR annotations)
                    continue
                parts = line.strip().split("\t")
                refseq = parts[0]
                refseq_endpts[refseq] = list(map(int, parts[3:]))  # Store the start/end points as integers
        return refseq_endpts

    def abs_distance(self, mrna_pos, endpts):
        """计算 mRNA 相对于每个区域的绝对距离"""
        return [mrna_pos - point for point in endpts]

    def rel_distance(self, mrna_pos, endpts):
        """计算 mRNA 相对于 5'UTR, CDS 和 3'UTR 的相对位置"""
        u5_st, u5_end, cds_st, cds_end, u3_st, u3_end = endpts
        if u5_st <= mrna_pos <= u5_end:
            return (mrna_pos - u5_st + 1) / (u5_end - u5_st + 1)
        elif cds_st <= mrna_pos <= cds_end:
            return (mrna_pos - cds_st + 1) / (cds_end - cds_st + 1) + 1
        elif u3_st <= mrna_pos <= u3_end:
            return (mrna_pos - u3_st + 1) / (u3_end - u3_st + 1) + 2
        return None  # if mrna_pos is outside of all regions

    def size_features(self, endpts):
        """计算 5'UTR, CDS 和 3'UTR 的大小"""
        u5_st, u5_end, cds_st, cds_end, u3_st, u3_end = endpts
        return [u5_end - u5_st, cds_end - cds_st, u3_end - u3_st]

    def process_bed_file(self):
        """处理 BED 文件，计算每个基因的相关信息"""
        excluded = 0
        total = 0

        with open(self.bed_file, 'r') as f , open(self.output, 'w') as out:
            out.write("chr\tcoord\tgene_name\trefseqID\trel_location\tutr5_st\tutr5_end\tcds_st\tcds_end\tutr3_st\tutr3_end\tutr5_size\tcds_size\tutr3_size\n")
            for line in f:
                if line.startswith("#"):  # Skip comment lines
                    continue
                total += 1
                parts = line.strip().split("\t")
                mrna_meta = parts[4]
                refseq, gname, region, mrna_pos = mrna_meta.split("|")

                if refseq not in self.refseq_endpts:
                    excluded += 1
                    continue

                endpts = self.refseq_endpts[refseq]
                mrna_pos = int(mrna_pos)

                abs_dist = self.abs_distance(mrna_pos, endpts)
                rel_dist = self.rel_distance(mrna_pos, endpts)
                feature_sizes = self.size_features(endpts)

                out.write(f"{parts[0]}\t{parts[2]}\t{gname}\t{refseq}\t{rel_dist}\t" + "\t".join(map(str, abs_dist)) + "\t" + "\t".join(map(str, feature_sizes)) + "\n")

        print(f"** Total bedfile lines: {total}")
        print(f"** Excluded (features not well defined or gene not found in gtf): {excluded}")
