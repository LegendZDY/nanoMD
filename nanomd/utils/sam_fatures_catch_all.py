#!/usr/bin/env python3
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20200603
Description:
This Script of ONTrnaDirectSeq from Legendzdy.  Of course, this means there's a possibility
for other ways.Use at your own discretion.
"""
import argparse
import re


def random_forest_data(samFile, base, fout):
    """
    This function will organize a dict of features.
    :param samFile: Pandas DataFrame to subset the two classes from
    :param fout: outfile
    """
    transcript_lengths = {}
    S = re.compile(r'\d+S')
    if base == "U":
        T = re.compile(r'T|t')
    elif base == "A":
        T = re.compile(r'A|a')
    elif base == "G":
        T = re.compile(r'G|g')
    elif base == "C":
        T = re.compile(r'C|c')
    else:
        print("Plase choose one of ATGC ")
    M = re.compile(r'\d+M')
    cstr = re.compile(r'cs:Z\S+')
    identical_block = re.compile(r':[0-9]+')
    deltion_block = re.compile(r'[\-][A-Za-z]+')
    insertion_block = re.compile(r'[\+][A-Za-z]+')
    substituted_ref_read = re.compile(r'\*[a-z][a-z]')

    with open(fout, "w") as fo:
        with open(samFile, "r") as fi:
            if base == "U":
                fo.write(
                    "read,UmisRate,UmisQrate,t_a_rate,t_a_Qrate,t_c_rate,t_c_Qrate,t_g_rate,t_g_Qrate,del_rate\n")
            elif base == "A":
                fo.write(
                    "read,AmisRate,AmisQrate,a_t_rate,a_t_Qrate,a_c_rate,a_c_Qrate,a_g_rate,a_g_Qrate,del_rate\n")
            elif base == "G":
                fo.write(
                    "read,GmisRate,GmisQrate,g_a_rate,g_a_Qrate,g_c_rate,g_c_Qrate,g_t_rate,g_t_Qrate,del_rate\n")
            elif base == "C":
                fo.write(
                    "read,CmisRate,CmisQrate,c_a_rate,c_a_Qrate,c_t_rate,c_t_Qrate,c_g_rate,c_g_Qrate,del_rate\n")
            else:
                print("Plase choose one of AUGC ")
            while True:
                try:
                    line = next(fi).strip()
                    if line.startswith('@'):
                        if line.startswith('@SQ'):
                            name = line[line.find('SN:') + 3:].split()[0]
                            length = float(
                                line[line.find('LN:') + 3:].split()[0])
                            transcript_lengths[name] = length
                        continue

                    if line.split()[1] == "0" or line.split()[1] == "16":
                        cs = line[cstr.search(line).span()[
                            0] + 5:cstr.search(line).span()[1]]
                        read, transcript, cigar, quality, sequence, Q = line.split()[0], line.split(
                        )[2], line.split()[5], int(line.split()[4]), line.split()[9], line.split()[10]
                        matchLen = sum([int(i[:-1]) for i in M.findall(cigar)])
                        coverage = matchLen/transcript_lengths[transcript]
                        if coverage < 0.8:
                            continue

                        readStart, UmisNum, UmisQ, delNum, t_a, t_c, t_g, t_aQ, t_cQ, t_gQ, matchseq = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ""
                        if S.match(cigar) != None:
                            numS = int(cigar[S.match(cigar).span()[
                                       0]:S.match(cigar).span()[1] - 1])
                            readStart = numS
                        while len(cs) != 0:
                            if identical_block.match(cs) != None:
                                idenum = int(cs[identical_block.match(cs).span()[
                                             0]+1:identical_block.match(cs).span()[1]])
                                if len(matchseq) == 0:
                                    matchseq += str(sequence[readStart:idenum])
                                else:
                                    matchseq += str(
                                        sequence[readStart:idenum+readStart])
                                readStart += idenum
                                cs = cs[identical_block.match(cs).span()[1]:]
                            elif deltion_block.match(cs) != None:
                                delStr = cs[deltion_block.match(cs).span(
                                )[0]+1:deltion_block.match(cs).span()[1]]
                                cs = cs[deltion_block.match(cs).span()[1]:]
                                if len(T.findall(delStr)) != 0:
                                    delNum += len(T.findall(delStr))
                            elif insertion_block.match(cs) != None:
                                readStart += len(cs[insertion_block.match(cs).span()[
                                                 0]+1:insertion_block.match(cs).span()[1]])
                                cs = cs[insertion_block.match(cs).span()[1]:]
                            elif substituted_ref_read.match(cs) != None:
                                refStr = cs[substituted_ref_read.match(cs).span(
                                )[0]+1:substituted_ref_read.match(cs).span()[1]][0]
                                readStr = cs[substituted_ref_read.match(cs).span(
                                )[0]+1:substituted_ref_read.match(cs).span()[1]][1]
                                cs = cs[substituted_ref_read.match(cs).span()[
                                    1]:]
                                if refStr == "t" and base == "U":
                                    accuracy = 1 - \
                                        0.1 ** ((ord(Q[readStart-1]) - 33) / 10)
                                    UmisNum += 1
                                    UmisQ += accuracy
                                    if readStr == "a":
                                        t_a += 1
                                        t_aQ += accuracy
                                    elif readStr == "c":
                                        t_c += 1
                                        t_cQ += accuracy
                                    else:
                                        t_g += 1
                                        t_gQ += accuracy
                                elif refStr == "a" and base == "A":
                                    accuracy = 1 - \
                                        0.1 ** ((ord(Q[readStart - 1]) - 33) / 10)
                                    UmisNum += 1
                                    UmisQ += accuracy
                                    if readStr == "t":
                                        t_a += 1
                                        t_aQ += accuracy
                                    elif readStr == "c":
                                        t_c += 1
                                        t_cQ += accuracy
                                    else:
                                        t_g += 1
                                        t_gQ += accuracy
                                elif refStr == "g" and base == "G":
                                    accuracy = 1 - \
                                        0.1 ** ((ord(Q[readStart - 1]) - 33) / 10)
                                    UmisNum += 1
                                    UmisQ += accuracy
                                    if readStr == "a":
                                        t_a += 1
                                        t_aQ += accuracy
                                    elif readStr == "c":
                                        t_c += 1
                                        t_cQ += accuracy
                                    else:
                                        t_g += 1
                                        t_gQ += accuracy
                                elif refStr == "c" and base == "C":
                                    accuracy = 1 - \
                                        0.1 ** ((ord(Q[readStart - 1]) - 33) / 10)
                                    UmisNum += 1
                                    UmisQ += accuracy
                                    if readStr == "a":
                                        t_a += 1
                                        t_aQ += accuracy
                                    elif readStr == "t":
                                        t_c += 1
                                        t_cQ += accuracy
                                    else:
                                        t_g += 1
                                        t_gQ += accuracy
                            else:
                                print("error in cs")
                                break
                        UmatchNum = len(T.findall(matchseq))
                        read_mean_Q = sum(
                            [1 - 0.1 ** ((ord(Q[i]) - 33) / 10) for i in range(len(Q))])/len(Q)
                        if UmisNum != 0:
                            UmisRate = UmisNum/(delNum+UmisNum+UmatchNum)
                            UmisQrate = (UmisQ/UmisNum)*read_mean_Q
                            t_a_rate = t_a / (delNum+UmisNum+UmatchNum)
                            t_c_rate = t_c / (delNum+UmisNum+UmatchNum)
                            t_g_rate = t_g / (delNum+UmisNum+UmatchNum)
                            if t_a != 0:
                                t_a_Qrate = (t_aQ/t_a)*read_mean_Q
                            else:
                                t_a_Qrate = 0
                            if t_c != 0:
                                t_c_Qrate = (t_cQ/t_c)*read_mean_Q
                            else:
                                t_c_Qrate = 0
                            if t_g != 0:
                                t_g_Qrate = (t_gQ/t_g)*read_mean_Q
                            else:
                                t_g_Qrate = 0
                        else:
                            UmisRate = 0
                            UmisQrate = 0
                            t_a_rate = 0
                            t_c_rate = 0
                            t_g_rate = 0
                            t_a_Qrate = 0
                            t_c_Qrate = 0
                            t_g_Qrate = 0
                        if delNum != 0:
                            del_rate = delNum/(delNum+UmisNum+UmatchNum)
                        else:
                            del_rate = 0

                        fo.write("{},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f}\n".format(
                            read, UmisRate, UmisQrate, t_a_rate, t_a_Qrate, t_c_rate, t_c_Qrate, t_g_rate, t_g_Qrate, del_rate))
                except StopIteration:
                    break


def main():
    parser = argparse.ArgumentParser(
        description='Extract length and quality information of sequences form fastq file')
    parser.add_argument('-s', '--input_sam', type=str, required=True,
                        help='sam file from minimap')
    parser.add_argument('-b', '--base', type=str, required=True,
                        help='choose one of AUGC')
    parser.add_argument('-o', '--out_file', type=str, required=False,
                        help='Path to an output file to be created,default= ./')
    args = parser.parse_args()
    random_forest_data(args.input_sam, args.base, args.out_file)


if __name__ == '__main__':
    main()
# random_forest_data("CTL2-0601.fastq.sam","test.csv")
