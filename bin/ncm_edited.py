#!/usr/bin/env python

### A modified version of the NGSCheckMate file ncm.py, 
### Original code by Sejoon Lee, Soo Lee & Alice E. Lee
### Modified by Simon Pearce, October 2022 and provided under the MIT licence, as is the original source code.
### This alters the R script by:
### * Making the dendrogram scale better as more samples are added
### * Ensuring a dendrogram is generated if only two samples are compared
### * Changing the background colour to white
### as well as incorporating two open pull requests
### https://github.com/parklab/NGSCheckMate/pull/31/commits/d1f5c809f61ac1439d7ae539a510da18fef5052e 
### by Gert Hulselmans
### and https://github.com/parklab/NGSCheckMate/pull/37/commits/7efa7e0ee4513b2ef5f2adb191ddc4a93b431c53
### by zmiimz.


import os
import math
import subprocess, time
import argparse
from argparse import RawTextHelpFormatter
from subprocess import call

global bed_file
global outdir
global outfilename
global temp_out
global testsamplename
global SAMTOOLS
global BCFTOOLS
global REF
global bam_list

glob_scores = dict()    #Whole score
feature_list = dict()   #Each Feature List
label = []              #Samples
features = []           #dbSNP features
mean_depth = dict()
real_depth = dict()
real_count = dict()
sum_file = dict()
out_tag = ""
pdf_tag = ""
Family_flag = False
Nonzero_flag = False


#Calculation of AVerages
def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

#Calulation of Pearson Correlation
def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)

## Need to be checked , n==0 case
    if n == 0 :
        return 0

    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    sqrt_xdiff2_ydiff2 = math.sqrt(xdiff2 * ydiff2)

    return diffprod / sqrt_xdiff2_ydiff2 if sqrt_xdiff2_ydiff2 != 0.0 else 0.0

# createDataSet
# base_dir : directory of files, bedFile: name of the bedFile
def createDataSetFromDir(base_dir, bedFile):
    for root, dirs, files in os.walk(base_dir):
        for file in files:
    	    if not file.endswith(".vcf"):
                continue

            link = root + '/' +  file
            f = open(link, "r")
            dbsnpf= open(bedFile,"r")
            depth = dict()
            depth[file] = 0
            real_count[file] = 0
            count = 0

            sum=dict()
            sum[file] = 0

            scores = dict()     # Scores of B-allel Frequencies
            #DBSNP ID collecting system
            for i in dbsnpf.readlines():
                temp = i.strip().split('\t')
                if temp[0].find("chr")!= -1:
                    ID = str(temp[0][3:]) + "_" + str(temp[2])
                else:   
                    ID = str(temp[0]) + "_" + str(temp[2])
                scores[ID] = 0
                count = count + 1
            
            ## 0618_samtools and haplotyper
            vcf_flag = 0

#            feature_list[file] = []
            score_set = dict()
            #VCF file PROCESSING  and Generation of features
            total = 0
            GVCF_samples = dict()
            for i in f.readlines():        
                if i.startswith("#"):
                    if i.find("DP4") != -1:
                        vcf_flag = 1
                    if i.find("#CHROM") != -1:
                        temp = i.strip().split('\t')
                        total=len(temp) - 9
                        if total != 1:
                            for sample_idx in range(0,total):
                                file = temp[sample_idx + 9]
                                GVCF_samples[temp[sample_idx + 9]] = []
                                score_set[temp[sample_idx + 9]] = dict()
                                depth[temp[sample_idx + 9]] = 0
                                real_count[temp[sample_idx + 9]] = 0
                                sum[temp[sample_idx + 9]] =0
                                feature_list[temp[sample_idx + 9]] = []
                        if total == 1:
                            feature_list[file] = []
                    continue

                temp = i.strip().split('\t')
              ## ID in BED file only 
                if temp[0].find("chr")!= -1:
                    ID = str(temp[0][3:]) + "_" + str(temp[1])
                else:   
                    ID = str(temp[0]) + "_" + str(temp[1])

                if ID not in scores:
                    continue

                if vcf_flag == 1:
                    values = temp[7].split(';')

                    if values[0].startswith("INDEL"):
                        continue

                    for j in values:
                        if j.startswith("DP4"):
                            readcounts = j.split(',')
                            readcounts[0] = readcounts[0][4:]
                            total_reads =(float(readcounts[0]) + float(readcounts[1]) + float(readcounts[2]) + float(readcounts[3])) 
                            score = 0
                            if total_reads > 0:
                                score = (float(readcounts[2]) + float(readcounts[3])) / total_reads
                                real_count[file] = real_count[file] + 1  

                            depth[file] = depth[file] + total_reads
                          
                            if ID in scores:
                                feature_list[file].append(ID)
                                scores[ID]= score
                                sum[file] = sum[file] + float(readcounts[2]) + float(readcounts[3])
                elif total == 1 and vcf_flag == 0:
                    format = temp[8].split(':')  ##Format
                    AD_idx = -1
                    DP_idx = -1
                    for idx in range(0,len(format)):
                        if format[idx] == "AD":
                            AD_idx = idx
                        elif format[idx] == "DP":
                            DP_idx = idx
                    if AD_idx == -1:
                        continue
                    if DP_idx == -1:
                        continue
                    idx = 9
                    values = temp[idx].split(":")
                    readcounts = values[AD_idx].split(',')

                    if float(readcounts[0]) + float(readcounts[1]) < 0.5:
                        score =0
                    else:
                        score = float(readcounts[1])/ (float(readcounts[0]) + float(readcounts[1]))
                    depth[file] = depth[file] + float(values[DP_idx])
                    if float(values[DP_idx]) > 0:
                        real_count[file] = real_count[file] + 1
                    
                    if ID in scores:
                        feature_list[file].append(ID)
                        scores[ID]= score  ##from here!
                        sum[file] = sum[file] + float(readcounts[1])
                else:  ###### Haplotyper or other VCF
                    format = temp[8].split(':')  ##Format
                    AD_idx = -1
                    DP_idx = -1
                    for idx in range(0,len(format)):
                        if format[idx] == "AD":
                            AD_idx = idx
                        elif format[idx] == "DP":
                            DP_idx = idx
                    if AD_idx == -1:
                        continue
                    if DP_idx == -1:
                        continue
                    idx = 9
                    for file in GVCF_samples:
                        values = temp[idx].split(":")
                        if len(values) < len(format):
                            score = 0
                            idx = idx + 1
                            continue

                        readcounts = values[AD_idx].split(',')

                        if float(readcounts[0]) + float(readcounts[1]) < 0.5:
                            score =0
                        else:
                            score = float(readcounts[1])/ (float(readcounts[0]) + float(readcounts[1]))
                        depth[file] = depth[file] + float(values[DP_idx])
                        if float(values[DP_idx]) > 0:
                            real_count[file] = real_count[file] + 1
                        if ID in scores:
                            feature_list[file].append(ID)
                            score_set[file][ID]= score   ##from here!
                            sum[file] = sum[file] + float(readcounts[1])

                        idx = idx + 1

## TOTAL is not 1 or total is 1 cases
            if total == 1:
                mean_depth[file] = depth[file] / float(count)
                real_depth[file] = depth[file] / float(real_count[file])
                sum_file[file] = sum[file]

                for key in features:
                    if glob_scores.has_key(file):
                        glob_scores[file].append(scores[key])
                    else:
                        glob_scores[file] = [scores[key]]    
            else:
                for file in GVCF_samples:
                    mean_depth[file] = depth[file] / float(count)
                    real_depth[file] = depth[file] / float(real_count[file])
                    sum_file[file] = sum[file]

                    for key in features:
                        if key not in score_set[file]:
                            score_set[file][key] = 0
                        if glob_scores.has_key(file):
                            glob_scores[file].append(score_set[file][key])
                        else:
                            glob_scores[file] = [score_set[file][key]]    
            dbsnpf.close()
            f.close()

    for key in sorted(glob_scores):
        label.append(key)

#create dataset from the VCF list files
def createDataSetFromList(base_list, bedFile):
    base_F = open(base_list,'r')
    for line in base_F.readlines():
        link = line.strip()
        f = open(link, "r")
        dbsnpf= open(bedFile,"r")
        file = os.path.basename(link)
        depth = dict()
        depth[file] = 0
        real_count[file] = 0
        count = 0

        sum=dict()
        sum[file] = 0

        scores = dict()     # Scores of B-allel Frequencies
        #DBSNP ID collecting system
        for i in dbsnpf.readlines():
            temp = i.strip().split('\t')
            if temp[0].find("chr")!= -1:
                ID = str(temp[0][3:]) + "_" + str(temp[2])
            else:   
                ID = str(temp[0]) + "_" + str(temp[2])
            scores[ID] = 0
            count = count + 1
        
        ## 0618_samtools and haplotyper
        vcf_flag = 0

#            feature_list[file] = []
        score_set = dict()
        #VCF file PROCESSING  and Generation of features
        total = 0
        GVCF_samples = dict()
        for i in f.readlines():        
            if i.startswith("#"):
                if i.find("DP4") != -1:
                    vcf_flag = 1
                if i.find("#CHROM") != -1:
                    temp = i.strip().split('\t')
                    total=len(temp) - 9
                    if total != 1:
                        for sample_idx in range(0,total):
                            file = temp[sample_idx + 9]
                            GVCF_samples[temp[sample_idx + 9]] = []
                            score_set[temp[sample_idx + 9]] = dict()
                            depth[temp[sample_idx + 9]] = 0
                            real_count[temp[sample_idx + 9]] = 0
                            sum[temp[sample_idx + 9]] =0
                            feature_list[temp[sample_idx + 9]] = []
                    if total == 1:
                        feature_list[file] = []
                continue

            temp = i.strip().split('\t')
          ## ID in BED file only 
            if temp[0].find("chr")!= -1:
                ID = str(temp[0][3:]) + "_" + str(temp[1])
            else:   
                ID = str(temp[0]) + "_" + str(temp[1])

            if ID not in scores:
                continue

            if vcf_flag == 1:
                values = temp[7].split(';')

                if values[0].startswith("INDEL"):
                    continue

                for j in values:
                    if j.startswith("DP4"):
                        readcounts = j.split(',')
                        readcounts[0] = readcounts[0][4:]
                        total_reads =(float(readcounts[0]) + float(readcounts[1]) + float(readcounts[2]) + float(readcounts[3])) 
                        score = 0
                        if total_reads > 0:
                            score = (float(readcounts[2]) + float(readcounts[3])) / total_reads
                            real_count[file] = real_count[file] + 1  

                        depth[file] = depth[file] + total_reads
                       
                        if ID in scores:
                            feature_list[file].append(ID)
                            scores[ID]= score
                            sum[file] = sum[file] + float(readcounts[2]) + float(readcounts[3])
            elif total == 1 and vcf_flag == 0:
                format = temp[8].split(':')  ##Format
                AD_idx = -1 
                DP_idx = -1
                for idx in range(0,len(format)):
                    if format[idx] == "AD":
                        AD_idx = idx
                    elif format[idx] == "DP":
                        DP_idx = idx
                if AD_idx == -1:
                    continue
                if DP_idx == -1:
                    continue
                idx = 9
                values = temp[idx].split(":")
                readcounts = values[AD_idx].split(',')

                if float(readcounts[0]) + float(readcounts[1]) < 0.5:
                    score =0
                else:
                    score = float(readcounts[1])/ (float(readcounts[0]) + float(readcounts[1]))
                depth[file] = depth[file] + float(values[DP_idx])
                if float(values[DP_idx]) > 0:
                    real_count[file] = real_count[file] + 1
                if ID in scores:
                    feature_list[file].append(ID)
                    scores[ID]= score  ##from here!
                    sum[file] = sum[file] + float(readcounts[1])
            else:  ###### Haplotyper or other VCF
                format = temp[8].split(':')  ##Format
                AD_idx = -1
                DP_idx = -1
                for idx in range(0,len(format)):
                    if format[idx] == "AD":
                        AD_idx = idx
                    elif format[idx] == "DP":
                        DP_idx = idx
                if AD_idx == -1:
                    continue
                if DP_idx == -1:
                    continue
                idx = 9
                for file in GVCF_samples:
                    values = temp[idx].split(":")
                    if len(values) < len(format):
                        score = 0
                        idx = idx + 1
                        continue

                    readcounts = values[AD_idx].split(',')

                    if float(readcounts[0]) + float(readcounts[1]) < 0.5:
                        score =0
                    else:
                        score = float(readcounts[1])/ (float(readcounts[0]) + float(readcounts[1]))
                    depth[file] = depth[file] + float(values[DP_idx])
                    if float(values[DP_idx]) > 0:
                        real_count[file] = real_count[file] + 1                   
 
                    if ID in scores:
                        feature_list[file].append(ID)
                        score_set[file][ID]= score   ##from here!
                        sum[file] = sum[file] + float(readcounts[1])

                    idx = idx + 1

## TOTAL is not 1 or total is 1 cases
        if total == 1:
            mean_depth[file] = depth[file] / float(count)
            real_depth[file] = depth[file] / float(real_count[file])
            sum_file[file] = sum[file]

            for key in features:
                if glob_scores.has_key(file):
                    glob_scores[file].append(scores[key])
                else:
                    glob_scores[file] = [scores[key]]    
        else:
            for file in GVCF_samples:
                mean_depth[file] = depth[file] / float(count)
                real_depth[file] = depth[file] / float(real_count[file])
                sum_file[file] = sum[file]

                for key in features:
                    if key not in score_set[file]:
                        score_set[file][key] = 0
                    if glob_scores.has_key(file):
                        glob_scores[file].append(score_set[file][key])
                    else:
                        glob_scores[file] = [score_set[file][key]]    
        dbsnpf.close()
        f.close()

    for key in sorted(glob_scores):
        label.append(key)


def createDataSetFromDir_TEST(base_dir, bedFile,order):
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if not file.endswith(".vcf"):
                continue

            link = root + '/' +  file
            f = open(link, "r")
            dbsnpf= open(bedFile,"r")
            depth = 0
            count = 0

            sum = 0

            scores = dict()     # Scores of B-allel Frequencies
            #DBSNP ID collecting system
            for i in dbsnpf.readlines():
                temp = i.strip().split('\t')
                ID = str(temp[0])+"_"+str(temp[2])
                scores[ID] = 0
                count = count + 1

            file = file + "_" + order
            feature_list[file] = []
            #VCF file PROCESSING  and Generation of features
            for i in f.readlines():
                if i.startswith("#"):
                    continue

                temp = i.split('\t')
                values = temp[7].split(';')

                if values[0].startswith("INDEL"):
                    continue

                for j in values:
                    if j.startswith("DP4"):
                        readcounts = j.split(',')
                        readcounts[0] = readcounts[0][4:]
                        score = (float(readcounts[2]) + float(readcounts[3])) / (float(readcounts[0]) + float(readcounts[1]) + float(readcounts[2]) + float(readcounts[3]))
                        depth = depth + (float(readcounts[0]) + float(readcounts[1]) + float(readcounts[2]) + float(readcounts[3]))

                        ID = str(temp[0]) + "_" + str(temp[1])
                        feature_list[file].append(ID)
                        scores[ID]= score
                        sum = sum + float(readcounts[2]) + float(readcounts[3])

            mean_depth[file] = depth / count
            sum_file[file] = sum

            for key in features:
                if glob_scores.has_key(file):
                    glob_scores[file].append(scores[key])
                else:
                    glob_scores[file] = [scores[key]]

            dbsnpf.close()
            f.close()

    for key in sorted(glob_scores):
        label.append(key)

#create dataset from the VCF list files
def createDataSetFromList_TEST(base_list, bedFile,order):
    base_F = open(base_list,'r')
    for line in base_F.readlines():
        link = line.strip()
        f = open(link, "r")
        dbsnpf= open(bedFile,"r")
        file = link[link.rindex("/")+1:]
        depth = 0
        count = 0

        sum = 0

        scores = dict()     # Scores of B-allel Frequencies
        #DBSNP ID collecting system
        for i in dbsnpf.readlines():
            temp = i.strip().split('\t')
            ID = str(temp[0])+"_"+str(temp[2])
            scores[ID] = 0
            count = count + 1

        file = file + "_" + order 
        feature_list[file] = []
        #VCF file PROCESSING  and Generation of features
        for i in f.readlines():
            if i.startswith("#"):
                continue

            temp = i.split('\t')
            values = temp[7].split(';')

            if values[0].startswith("INDEL"):
                continue

            for j in values:
                if j.startswith("DP4"):
                    readcounts = j.split(',')
                    readcounts[0] = readcounts[0][4:]
                    score = (float(readcounts[2]) + float(readcounts[3])) / (float(readcounts[0]) + float(readcounts[1]) + float(readcounts[2]) + float(readcounts[3]))
                    depth = depth + (float(readcounts[0]) + float(readcounts[1]) + float(readcounts[2]) + float(readcounts[3]))

                    ID = str(temp[0]) + "_" + str(temp[1])
                    feature_list[file].append(ID)
                    scores[ID]= score
                    sum = sum + float(readcounts[2]) + float(readcounts[3])

        mean_depth[file] = depth / count
        sum_file[file] = sum

        for key in features:
            if glob_scores.has_key(file):
                glob_scores[file].append(scores[key])
            else:
                glob_scores[file] = [scores[key]]

        dbsnpf.close()
        f.close()

    for key in sorted(glob_scores):
        label.append(key)


# kNN based classification
def clustering(K):
    altFreqList = []
    keyList = []
    Pos_count = 0

    for key in sorted(glob_scores):
        altFreqList.append(glob_scores[key])
        keyList.append(key)

    dataSetSize = len(altFreqList)

    sum = 0
    othersum = 0
    for target in range(0,dataSetSize):
        dist = []
        pheno = []
        # comparison to the other samples based on BASE sample
        base = altFreqList[target]
        tempA = set(feature_list[keyList[target]])

        # calculate eucladian distance between two samples
        for i in range(0, dataSetSize):

#            IsdiffPhenotype = 0.0
            comparison = altFreqList[i]

            tempB = set(feature_list[keyList[i]])

            selected_feature = tempA.intersection(tempB)

#            IsdiffPhenotype = (2*len(selected_feature))/(len(tempA) + len(tempB))

            vecA = []
            vecB = []

            idx = 0
            for k in features:
                if k in selected_feature:
                    vecA.append(base[idx])
                    vecB.append(comparison[idx])
                idx = idx + 1

            distance = pearson_def(vecA, vecB)
            dist.append(distance)
#            pheno.append(IsdiffPhenotype)

        orderCount = 0
        while (orderCount < K):
            max_value = sorted(dist)[-2-orderCount]
            max_indice = dist.index(max_value)
            sum = sum + max_value
            Pos_count =  Pos_count + 1
            outPOS=str(label[target]) +  "\tmatched to\t" + str(label[max_indice])+ "\tscore=\t" + str(max_value)
            print outPOS
            #POS_F.write(outPOS + "\n")
            orderCount = orderCount + 1

#    print sum/Pos_count

#OLD version
def classify(T):
    altFreqList = []
    keyList = []
    Pos_count = 0
    Neg_count = 0

    POS_F = open("/data/users/sjlee/valid_qc/WGS/SNP/results/TEST2_POS_SNP.txt",'w')
    NEG_F = open("/data/users/sjlee/valid_qc/WGS/SNP/results/TEST2_NEG_SNP.txt",'w')

    for key in sorted(glob_scores):
        altFreqList.append(glob_scores[key])
        keyList.append(key)

    dataSetSize = len(altFreqList)

    sum = 0
    othersum = 0
    for target in range(0,dataSetSize):
        dist = []
        pheno = []
        # comparison to the other samples based on BASE sample
        base = altFreqList[target]
        tempA = set(feature_list[keyList[target]])

        # calculate eucladian distance between two samples
        for i in range(0, dataSetSize):

            IsdiffPhenotype = 0.0
            comparison = altFreqList[i]

            tempB = set(feature_list[keyList[i]])

            selected_feature = tempA.intersection(tempB)

            IsdiffPhenotype = (2*len(selected_feature))/(len(tempA) + len(tempB))

            vecA = []
            vecB = []

            idx = 0
            for k in features:
                if k in selected_feature:
                    vecA.append(base[idx])
                    vecB.append(comparison[idx])
                idx = idx + 1

            distance = pearson_def(vecA, vecB)
            dist.append(distance)
            pheno.append(IsdiffPhenotype)

        for value in sorted(dist)[0:-2]:
            if abs((Tmean-value)/Tstd) < abs((Fmean-value)/Fstd):
                max_value = value
                max_indice = dist.index(max_value)
                td = array(dist)
                sum = sum + max_value
                Pos_count =  Pos_count + 1
                outPOS=str(label[target]) +  "\tmatched to\t" + str(label[max_indice])+ "\tscore=\t" + str(max_value) + "\tdiff=\t" + str(pheno[max_indice])
                POS_F.write(outPOS + "\n")
            else:
                max_value = value
                max_indice = dist.index(max_value)
                othersum = othersum + max_value
                Neg_count = Neg_count + 1
                outNEG=str(label[target]) +  "\tmatched to\t" + str(label[max_indice])+ "\tscore=\t" + str(max_value) + "\tdiff=\t" + str(pheno[max_indice])
                NEG_F.write(outNEG + "\n")


    print sum/Pos_count
    print othersum/Neg_count

    POS_F.close()
    NEG_F.close()


def classifyNV(vec2Classify, p0Vec, p0S, p1Vec, p1S):
    if abs(p0Vec - vec2Classify) - p0S > abs(p1Vec - vec2Classify) - p1S:
        return abs((abs(p0Vec - vec2Classify) - p0S )/ (abs(p1Vec - vec2Classify) -  p1S )), 1
    else: 
        return abs((abs(p0Vec - vec2Classify) - p0S) / (abs(p1Vec - vec2Classify)  -  p1S)), 0  

#    if depth < 5:
#        if (vec2Classify >= (p1Vec - p1S)):
#            return (abs(p0Vec - vec2Classify) / p0S )/ (abs(p1Vec - vec2Classify)/ p1S ), 1
#        else:
#            return (abs(p0Vec - vec2Classify) / p0S) / (abs(p1Vec - vec2Classify)/ p1S), 0
#    else:
#        if (abs(p0Vec - vec2Classify) / p0S > abs(p1Vec - vec2Classify)/ p1S):
#            return (abs(p0Vec - vec2Classify) / p0S )/ (abs(p1Vec - vec2Classify)/ p1S ), 1
#        else:
#            return (abs(p0Vec - vec2Classify) / p0S) / (abs(p1Vec - vec2Classify)/ p1S), 0


def trainNV(trainMatrix,trainCategory):
    numTrainDocs = len(trainMatrix)    # #of traning samples

    p1List = []
    p0List = []

    for i in range(numTrainDocs):
        if trainCategory[i] == 1:
            p1List.append(trainMatrix[i])
        else:
            p0List.append(trainMatrix[i])

    return mean(p1List),std(p1List), mean(p0List),std(p0List)


def calAUC(predStrengths, classLabels):
    ySum = 0.0 #variable to calculate AUC
    cur = (1.0,1.0) #cursor
    numPosClas = sum(array(classLabels)==1.0)
    yStep = 1/float(numPosClas); xStep = 1/float(len(classLabels)-numPosClas)
    sortedIndicies = predStrengths.argsort()#get sorted index, it's reverse
    #loop through all the values, drawing a line segment at each point
    for index in sortedIndicies.tolist()[0]:
        if classLabels[index] == 1:
            delX = 0; delY = yStep;
        else:
            delX = xStep; delY = 0;
            ySum += cur[1]
        cur = (cur[0]-delX,cur[1]-delY)
    return ySum*xStep

#def plotROC(predStrengths, classLabels):
#    import matplotlib.pyplot as plt
#    cur = (1.0,1.0) #cursor
#    ySum = 0.0 #variable to calculate AUC
#    numPosClas = sum(array(classLabels)==1.0)
#    yStep = 1/float(numPosClas); xStep = 1/float(len(classLabels)-numPosClas)
#    sortedIndicies = predStrengths.argsort()#get sorted index, it's reverse
#    fig = plt.figure()
#    fig.clf()
#    ax = plt.subplot(111)
#    #loop through all the values, drawing a line segment at each point
#    for index in sortedIndicies.tolist()[0]:
#        if classLabels[index] == 1:
#            delX = 0; delY = yStep;
#        else:
#            delX = xStep; delY = 0;
#            ySum += cur[1]
#        #draw line from cur to (cur[0]-delX,cur[1]-delY)
#        ax.plot([cur[0],cur[0]-delX],[cur[1],cur[1]-delY], c='b')
#        cur = (cur[0]-delX,cur[1]-delY)
#    ax.plot([0,1],[0,1],'b--')
#    plt.xlabel('False positive rate'); plt.ylabel('True positive rate')
#    plt.title('ROC curves')
#    ax.axis([0,1,0,1])
#    plt.show()
#    print "the Area Under the Curve is: ",ySum*xStep

def getPredefinedModel(depth):
     if Family_flag:
         if depth > 10:
             return 0.874611,0.022596,0.644481,0.020908
         elif depth > 5:
             return 0.785312,0.021318,0.596133,0.022502
         elif depth > 2:
             return 0.650299,0.019252,0.5346,0.020694
         elif depth > 1:
             return 0.578582,0.018379,0.495017,0.021652
         elif depth > 0.5:
             return 0.524757,0.023218,0.465653,0.027378
         else:
    #         print "Warning: Sample region depth is too low < 1"
             return 0.524757,0.023218, 0.465653, 0.027378
     else:
         if depth > 10:
             return 0.874546, 0.022211, 0.310549, 0.060058
         elif depth > 5:
             return 0.785249,0.021017, 0.279778, 0.054104
         elif depth > 2:
             return 0.650573, 0.018699,0.238972, 0.047196
         elif depth > 1:
             return 0.578386,0.018526, 0.222322, 0.041186
         elif depth > 0.5:
             return 0.529327,0.025785, 0.217839, 0.040334
         else:
    #         print "Warning: Sample region depth is too low < 1"
             return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 30:
#         return 0.874546, 0.022211, 0.310549, 0.060058
#     elif depth > 10:
#         return 0.785249,0.021017, 0.279778, 0.054104
#     elif depth > 5:
#         return 0.650573, 0.018699,0.238972, 0.047196
#     elif depth > 2:
#         return 0.578386,0.018526, 0.222322, 0.041186
#     elif depth > 1:
#         return 0.529327,0.025785, 0.217839, 0.040334
#     else:
#         print "Warning: Sample region depth is too low < 1"
#         return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 0.1:
#        return 0.0351* depth + 0.5538, 0.02, 0.009977*depth + 0.216978, 0.045
#     else:
#        print "too low depth"
#        return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 0.5:
#        return 0.06315* (math.log(depth)) + 0.64903, 0.046154, 0.0005007*depth + 0.3311504,0.12216
#     else:
#        return 0.62036, 0.046154, 0.31785, 0.12216

def getPredefinedModel_F(depth):
     if depth > 10:
         return 0.874546, 0.022211, 0.620549, 0.060058
     elif depth > 5:
         return 0.785249,0.021017, 0.609778, 0.054104
     elif depth > 2:
         return 0.650573, 0.018699,0.548972, 0.047196
     elif depth > 1:
         return 0.578386,0.018526, 0.502322, 0.041186
     elif depth > 0.5:
         return 0.529327,0.025785, 0.457839, 0.040334
     else:
#         print "Warning: Sample region depth is too low < 1"
         return 0.529327,0.025785, 0.457839, 0.040334
#     if depth > 30:
#         return 0.874546, 0.022211, 0.310549, 0.060058
#     elif depth > 10:
#         return 0.785249,0.021017, 0.279778, 0.054104
#     elif depth > 5:
#         return 0.650573, 0.018699,0.238972, 0.047196
#     elif depth > 2:
#         return 0.578386,0.018526, 0.222322, 0.041186
#     elif depth > 1:
#         return 0.529327,0.025785, 0.217839, 0.040334
#     else:
#         print "Warning: Sample region depth is too low < 1"
#         return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 0.1:
#        return 0.0351* depth + 0.5538, 0.02, 0.009977*depth + 0.216978, 0.045
#     else:
#        print "too low depth"
#        return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 0.5:
#        return 0.06315* (math.log(depth)) + 0.64903, 0.046154, 0.0005007*depth + 0.3311504,0.12216
#     else:
#        return 0.62036, 0.046154, 0.31785, 0.12216


def classifying():
    AUCs =[]

    wholeFeatures = 50

    temp =[]

    altFreqList = []
    keyList = []

    for key in sorted(glob_scores):
        altFreqList.append(glob_scores[key])
        keyList.append(key)

    dataSetSize = len(altFreqList)

    filter_list = []

    for i in range(0, dataSetSize):
        for j in range(0, dataSetSize):
            if i!=j:
                if keyList[j] not in filter_list:
                    temp.append([keyList[i],keyList[j]])
        filter_list.append(keyList[i])

    for iterations in range(49,wholeFeatures):

        samples = []
        numFeatures = iterations

        count = 0

        for i in range(0,len(temp)):
            tempA = set(feature_list[temp[i][0].strip()])
            tempB = set(feature_list[temp[i][1].strip()])

            selected_feature = tempA.intersection(tempB)

            vecA = []
            vecB = []

            idx = 0
            for k in features:
                if k in selected_feature:
                    vecA.append(glob_scores[temp[i][0].strip()][idx])
                    vecB.append(glob_scores[temp[i][1].strip()][idx])
                idx = idx + 1

            distance = pearson_def(vecA, vecB)
            samples.append(distance)

        predStrength = []
        training_flag =0
    ####0715 Append

        output_matrix_f = open(outdir + "/" + out_tag + "_output_corr_matrix.txt","w")
        output_matrix = dict()
        
        if out_tag!="stdout":
        	out_f = open(outdir + "/" + out_tag + "_all.txt","w")
        	out_matched = open(outdir + "/" + out_tag + "_matched.txt","w")

        for i in range(0, len(keyList)):
            output_matrix[keyList[i]] = dict()
            for j in range(0,len(keyList)):
                output_matrix[keyList[i]][keyList[j]] = 0

        if training_flag == 1:
            #make training set
            for i in range(0,len(samples)):
                trainMatrix= []
                trainCategory = []
                for j in range(0, len(samples)):
                    if i==j:
                        continue
                    else:
                        trainMatrix.append(samples[j])
                        trainCategory.append(classLabel[j])
                #training samples in temp
                #p0V, p1V, pAb = trainNB0(array(trainMatrix),array(trainCategory))
                p1V,p1S, p0V, p0S = trainNV(array(trainMatrix),array(trainCategory))
                result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
                if result[1] == 1:
                    print str(temp[i][0]) + '\tsample is matched to\t',str(temp[i][1]),'\t', samples[i]
                predStrength.append(result[0])
        else :
            for i in range(0,len(samples)):
                depth = 0
                if Nonzero_flag: 
                    depth = min(real_depth[temp[i][0].strip()],real_depth[temp[i][1].strip()])
                else:
                    depth = min(mean_depth[temp[i][0].strip()],mean_depth[temp[i][1].strip()])
               
                p1V,p1S, p0V, p0S = getPredefinedModel(depth)
                result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
                if result[1] ==1:
                    output_matrix[temp[i][0].strip()][temp[i][1].strip()] = samples[i]
                    output_matrix[temp[i][1].strip()][temp[i][0].strip()] = samples[i]
                    if out_tag=="stdout":
                        print str(temp[i][0]) + '\tmatched\t',str(temp[i][1]),'\t', round(samples[i],4),'\t',round(depth,2)
                    else :
                        out_f.write(str(temp[i][0]) + '\tmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                        out_matched.write(str(temp[i][0]) + '\tmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                else:
                    if out_tag=="stdout":
                        print str(temp[i][0]) + '\tunmatched\t',str(temp[i][1]),'\t', round(samples[i],4),'\t',round(depth,2)
                    else :
                        out_f.write(str(temp[i][0]) + '\tunmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                predStrength.append(result[0])
            #testing sample is samples
        output_matrix_f.write("sample_ID")
        for key in output_matrix.keys():
            if key.find(".vcf") != -1:
                output_matrix_f.write("\t" + key[0:key.index('.vcf')])
            else:
                output_matrix_f.write("\t" + key)
        output_matrix_f.write("\n")

#        for key in output_matrix.keys():
#            for otherkey in output_matrix[key].keys():
#                if output_matrix[key][otherkey] != 0:
#                    output_matrix[otherkey][key] = output_matrix[key][otherkey] 

        for key in output_matrix.keys():
            if key.find(".vcf") != -1:
                output_matrix_f.write(key[0:key.index('.vcf')])
            else:
                output_matrix_f.write(key)
            for otherkey in output_matrix.keys():
                output_matrix_f.write("\t" + str(output_matrix[key][otherkey]))
            output_matrix_f.write("\n")   
            
        output_matrix_f.close()         
        if out_tag!="stdout":
        	out_f.close()
        	out_matched.close()   



def classifying_test():
    AUCs =[]

    wholeFeatures = 50

    temp = []

    keyF = open(testsamplename,'r')
    temp =[]

    for k in keyF.readlines():
        keyfile = k.split(":")
        keyfile[0] = keyfile[0].strip() + "_1"
        keyfile[1] = keyfile[1].strip() + "_2"
        temp.append(keyfile)
    keyF.close()

    for iterations in range(49,wholeFeatures):

        samples = []
        numFeatures = iterations

        count = 0

        for i in range(0,len(temp)):
            tempA = set(feature_list[temp[i][0].strip()])
            tempB = set(feature_list[temp[i][1].strip()])

            selected_feature = tempA.intersection(tempB)
            
            vecA = []
            vecB = []
            
            idx = 0
            for k in features:
                if k in selected_feature:
                    vecA.append(glob_scores[temp[i][0].strip()][idx])
                    vecB.append(glob_scores[temp[i][1].strip()][idx])
                idx = idx + 1
            
            distance = pearson_def(vecA, vecB)
            samples.append(distance)
            
        predStrength = []
        training_flag =0
    ####0715 Append

        output_matrix_f = open(outdir + "/" + out_tag + "_output_corr_matrix.txt","w")
        output_matrix = dict()
        
        if out_tag!="stdout":
            out_f = open(outdir + "/" + out_tag + "_all.txt","w")
            out_matched = open(outdir + "/" + out_tag + "_matched.txt","w")

        for i in range(0, len(keyList)):
            output_matrix[keyList[i]] = dict()
            for j in range(0,len(keyList)):
                output_matrix[keyList[i]][keyList[j]] = 0

        if training_flag == 1:
            #make training set
            for i in range(0,len(samples)):
                trainMatrix= []
                trainCategory = []
                for j in range(0, len(samples)):
                    if i==j:
                        continue
                    else:
                        trainMatrix.append(samples[j])
                        trainCategory.append(classLabel[j])
                #training samples in temp
                #p0V, p1V, pAb = trainNB0(array(trainMatrix),array(trainCategory))
                p1V,p1S, p0V, p0S = trainNV(array(trainMatrix),array(trainCategory))
                result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
                if result[1] == 1:
                    print str(temp[i][0]) + '\tsample is matched to\t',str(temp[i][1]),'\t', samples[i]
                predStrength.append(result[0])
        else :
            for i in range(0,len(samples)):
                depth = min(mean_depth[temp[i][0].strip()],mean_depth[temp[i][1].strip()])
                p1V,p1S, p0V, p0S = getPredefinedModel(depth)
                result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
                if result[1] ==1:
                    output_matrix[temp[i][0].strip()][temp[i][1].strip()] = samples[i]
                    output_matrix[temp[i][1].strip()][temp[i][0].strip()] = samples[i]
                    if out_tag=="stdout":
                        print str(temp[i][0]) + '\tmatched\t',str(temp[i][1]),'\t', round(samples[i],4),'\t',round(depth,2)
                    else :
                        out_f.write(str(temp[i][0]) + '\tmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                        out_matched.write(str(temp[i][0]) + '\tmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                else:
                    if out_tag=="stdout":
                        print str(temp[i][0]) + '\tunmatched\t',str(temp[i][1]),'\t', round(samples[i],4),'\t',round(depth,2)
                    else :
                        out_f.write(str(temp[i][0]) + '\tunmatched\t' + str(temp[i][1])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                predStrength.append(result[0])
            #testing sample is samples
        output_matrix_f.write("sample_ID")
        for key in output_matrix.keys():
            output_matrix_f.write("\t" + key[0:key.index('.')])
        output_matrix_f.write("\n")

#        for key in output_matrix.keys():
#            for otherkey in output_matrix[key].keys():
#                if output_matrix[key][otherkey] != 0:
#                    output_matrix[otherkey][key] = output_matrix[key][otherkey] 

        for key in output_matrix.keys():
            output_matrix_f.write(key[0:key.index('.')])
            for otherkey in output_matrix.keys():
                output_matrix_f.write("\t" + str(output_matrix[key][otherkey]))
            output_matrix_f.write("\n")   
            
        output_matrix_f.close()         
        if out_tag!="stdout":
            out_f.close()
            out_matched.close()



def generate_R_scripts():
    r_file = open(outdir + "/r_script.r","w")
    if len(feature_list)==0:
       r_file.close()
    else :
       cmd = "output_corr_matrix <- read.delim(\"" + outdir +  "/" + out_tag + "_output_corr_matrix.txt\")\n"
       cmd = cmd + "data = output_corr_matrix\n"
       cmd = cmd + "d3 <- as.dist((1 - data[,-1]))\n"
       cmd = cmd + "clust3 <- hclust(d3, method = \"average\")\n"
       if len(feature_list) < 5:
           cmd = cmd + "pdf(\"" +outdir+ "/" + pdf_tag + ".pdf\", width=10, height=7)\n"
       else:
           cmd = cmd + "pdf(\"" +outdir+ "/" + pdf_tag + ".pdf\", width="+str(math.log10(7*len(feature_list))*10) +", height=7)\n"
       cmd = cmd + "op = par(bg = \"white\")\n"
       cmd = cmd + "par(plt=c(0.05, 0.95, 0.25, 0.9))\n"
       if len(feature_list) < 3:
           cmd = cmd + "plot(as.dendrogram(clust3), lwd = 2, lty = 1,cex=0.8, xlab=\"Samples\", sub = \"\",  ylab=\"Distance (1-Pearson correlation)\", axes = FALSE)\n"
       else:
           cmd = cmd + "plot(clust3, lwd = 2, lty = 1,cex=0.8, xlab=\"Samples\", sub = \"\",  ylab=\"Distance (1-Pearson correlation)\",hang = -1, axes = FALSE)\n"
       cmd = cmd + "axis(side = 2, at = seq(0, 1, 0.2), labels = FALSE, lwd = 2)\n"
       cmd = cmd + "mtext(seq(0, 1, 0.2), side = 2, at = seq(0, 1, 0.2), line = 1,   las = 2)\n"
       cmd = cmd + "dev.off()\n"
       r_file.write(cmd)
       r_file.close()



def run_R_scripts():
    command = "R CMD BATCH " + outdir + "/r_script.r"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()


def remove_internal_files():
    if outdir.find("*"):
        sys.exit()


    command = "rm -rf " + outdir + "/" + out_tag + "_output_corr_matrix.txt"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()
    command = "rm -rf " + outdir + "/r_script.r"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()
    command = "rm -rf " + outdir + "/r_script.r.Rout"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()


def getCallResult(command):
    fd_popen = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    (stdoutdata,stderrdata) = fd_popen.communicate()
    return stdoutdata,stderrdata


def run_mpileup():
    

    SAMTOOLS="samtools"
    BCFTOOLS="bcftools"
    REF=os.environ.get("NCM_REF",False)
    
    version =""
##version of samtools
    samtools_version = getCallResult(SAMTOOLS)
    for samtool_line in samtools_version:
        if samtool_line.find("Version") != -1:
            version_flag = 1
            version_line = samtool_line.split("\n")
            for version_tag in version_line:
                if version_tag.find("Version") != -1:
                    version_list = version_tag.split(" ")
                    version = version_list[1]
    print version

    for sample in bam_list:
        filename = sample.split("/")
        tag = filename[-1][0:filename[-1].rindex(".")]
        if version.startswith("0"):
            command = SAMTOOLS + " mpileup -I -uf " + REF + " -l " + bedFile + " " + sample + " | "  + BCFTOOLS + " view -cg - > " + outdir + "/" + tag  + ".vcf"
        else:
            command = SAMTOOLS + " mpileup -uf " + REF + " -l " + bedFile + " " + sample + " | "  + BCFTOOLS + " call -c > " + outdir + "/" + tag  + ".vcf"
        print command
        call(command,shell=True)
 #       proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
#        return_code = proc.wait()

def find_bam_list():
    for root, dirs, files in os.walk(base_dir):
        for file in files:
             if not file.endswith(".bam"):
                continue
             bam_list.append(root + "/" + file)

def get_bam_list():
    with open(base_list,'r') as F:
        for line in F.readlines():
             bam_list.append(line.strip())


def output_filter():
    success_set_M = []
    success_set_U = []
    failure_set_M = []
    failure_set_U = []

    with open(outdir + "/" + out_tag + "_all.txt","r") as F:
        for line in F.readlines():
            temp = line.strip().split('\t')
            
            sample1 = temp[0]
            sample2 = temp[2]
            
            match = temp[1]
            
            if match == "matched":
                if sample1[sample1.index("TCGA"):sample1.index("TCGA")+12] == sample2[sample2.index("TCGA"):sample2.index("TCGA")+12] :
                    success_set_M.append(line)
                else:
                    failure_set_M.append(line)
            elif match == "unmatched":
                if sample1[sample1.index("TCGA"):sample1.index("TCGA")+12] == sample2[sample2.index("TCGA"):sample2.index("TCGA")+12] :
                    failure_set_U.append(line)
                else:
                    success_set_U.append(line)        
              
    Matched_file = open(outdir + "/" + out_tag + "_matched.txt",'w') 

    for i in success_set_M:
        Matched_file.write(i)
    for i in failure_set_M:
        Matched_file.write(i)  
    
    Matched_file.close()

    problem_file = open(outdir + "/" + out_tag + "_problematic.txt",'w')

    for i in failure_set_M:
        problem_file.write(i)
    for i in failure_set_U:
        problem_file.write(i)

    problem_file.close()

    Summary_file = open(outdir + "/" + out_tag + "_summary.txt",'w')
    
 

    ## paired cluster - only failed things
    Summary_file.write("###########################################\n")
    Summary_file.write("###  Problematic clusters of same orgins ##\n")
    Summary_file.write("###########################################\n\n")

    cluster = dict()

    result_set = failure_set_M + success_set_M

    for line in result_set:
        temp = line.strip().split('\t')
        flag = 0
        for key in cluster:
            if temp[0] in cluster[key]:
                cluster[key].add(temp[2])
                flag = 1
                break
            elif temp[2] in cluster[key]:
                cluster[key].add(temp[0])
                flag = 1
                break
        
        if flag == 0:
            cluster[temp[0]] = set()
            cluster[temp[0]].add(temp[0])
            cluster[temp[0]].add(temp[2])
            
            
    count = 0 
    for key in cluster:
        temp_list = []
        flag = 0
        for data in cluster[key]:
            temp_list.append(data)
            sample1 = temp_list[0]
            ID = sample1[sample1.index("TCGA"):sample1.index("TCGA")+12]
            
            for sample1 in cluster[key]:
                if ID != sample1[sample1.index("TCGA"):sample1.index("TCGA")+12]:
                    flag = 1

              

        if flag == 1:
            count = count + 1
            Summary_file.write("Cluster " + str(count) + "\n")
              
            for data in cluster[key]:
                Summary_file.write(data + "\n")
            Summary_file.write("\n")

                
    ## Singleton
    Summary_file.write("\n")
    Summary_file.write("###########################################\n")
    Summary_file.write("############### Singleton #################\n")
    Summary_file.write("###########################################\n\n")

    final_set = set()
    filter_set = set()

    result_set = failure_set_U

    for line in result_set:
        temp = line.strip().split('\t')
        
        final_set.add(temp[0])
        final_set.add(temp[2])
        
        flag = 0
        for key in cluster:
            if temp[0] in cluster[key]:
                filter_set.add(temp[0])
            elif temp[2] in cluster[key]:
                filter_set.add(temp[2])
                


    for i in final_set.difference(filter_set):
        Summary_file.write(i + "\n")

    Summary_file.close()


if __name__ == '__main__':
    testsamplename = ""

    help = """
    Ensuring Sample Identity v1.0_modified
    Usage:   NGSCheckmate
    Desc.:   Input = the absolute path list of vcf files (samtools mpileup and bcftools)
             Output = Matched samples List
    Eg.  :   ncm.py -V -l /data/vcf_list.txt -bed /data/SNP_hg19.bed -O /data/output/ -N Matched_list
             ncm.py -V -d /data/vcf/ -bed /data/SNP_hg19.bed -O /data/output -N Matched_list
             ncm.py -B -d /data/bam/ -bed /data/SNP_hg19.bed -O /data/output -N Matched_list
             ncm.py -B -l /data/bam_list.txt -bed /data/SNP_hg19.bed -O /data/output/ -N Matched_list
    Sejoon Lee, Soo Lee, Eunjung Lee, 2015
    """

    parser = argparse.ArgumentParser(description=help, formatter_class=RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    datatype = parser.add_mutually_exclusive_group(required=True)
    datatype.add_argument('-B','--BAM',dest='BAM_type',action='store_true', help='BAM files to do NGS Checkmate')
    datatype.add_argument('-V','--VCF',dest='VCF_type',action='store_true', help='VCF files to do NGS Checkmate')
    parser.add_argument('-f','--family_cutoff',dest='family_cutoff',action='store_true', help='apply strict correlation threshold to remove family cases')

    group.add_argument('-l','--list',metavar='data_list',dest='datalist',action='store', help='data list')
    group.add_argument('-d','--dir',metavar='data_dir',dest='datadir',action='store', help='data directory')

    parser.add_argument('-bed','--bedfile',metavar='BED file',required=True,dest='bed_file',action='store', help='A bed file containing SNP polymorphisms')

    parser.add_argument('-t','--testsamplename',metavar='test_samplename',dest='testsamplename',action='store',help='file including test sample namses  with ":" delimeter (default : all combinations of samples), -t filename')
    parser.add_argument('-O','--outdir',metavar='output_dir',dest='outdir',action='store', help='directory name for temp and output files')
    parser.add_argument('-N','--outfilename',metavar='output_filename',dest='outfilename',action='store',default="output",help='OutputFileName ( default : output ), -N filename')
    parser.add_argument('-nz','--nonzero',dest='nonzero_read',action='store_true',help='Use non-zero mean depth of target loci as reference correlation. (default: Use mean depth of all target loci)')

#    parser.add_argument('-m','--method',metavar='METHOD',required=True,dest='method',action='store', choices={'clustering','classifying'},default='classifying', help='Method (Clustering, Classifying)')
#    parser.add_argument('-k','--knn',metavar='1,2,3,..',dest='KNN',action='store', default="1", help='K value for K-nearest neighbors clustering')
 #   parser.add_argument('-p','--pdf',metavar='PDFfile',dest='PDF_flag',action='store',default="otuput",help='output Prgramming R based clustering vector image of PDF file, -p filename')

    args=parser.parse_args()

    base_list = ""
    base_dir = ""

    if args.datalist != None :
        base_list = args.datalist
    if args.datadir != None :
        base_dir = args.datadir
        
    if args.family_cutoff:
        Family_flag=True
    if args.nonzero_read:
        Nonzero_flag=True

    outdir = args.outdir
    bedFile = args.bed_file
    out_tag = args.outfilename

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    key_order = open(bedFile,'r')
    for line in key_order.readlines():
        if line.startswith("#"):
            continue
        temp = line.strip().split('\t')
        if temp[0].find("chr")!= -1:
            features.append(str(temp[0][3:])+"_"+ str(temp[2]))
        else:   
            features.append(str(temp[0])+"_"+ str(temp[2]))
        

    if args.BAM_type != False :
        if args.datadir != None :
            bam_list = []
            find_bam_list()
        elif args.datalist != None :
            bam_list = []
            get_bam_list()
        run_mpileup()
        base_dir = outdir
        print "Generate Data Set from " + base_dir + "\nusing this bed file : " + bedFile
        if args.testsamplename != None:
            testsamplename = args.testsamplename
            createDataSetFromDir_TEST(base_dir,bedFile,"1")
            createDataSetFromDir_TEST(base_dir,bedFile,"2")
            classifying_test()
        else:
            createDataSetFromDir(base_dir,bedFile)
            classifying()
    elif args.VCF_type != False :
        if args.datadir != None :
            print "Generate Data Set from " + base_dir + "\nusing this bed file : " + bedFile
            if args.testsamplename != None:
                testsamplename = args.testsamplename
                createDataSetFromDir_TEST(base_dir,bedFile,"1")
                createDataSetFromDir_TEST(base_dir,bedFile,"2")
                classifying_test()
            else:
                createDataSetFromDir(base_dir,bedFile)
                classifying()
        elif args.datalist != None :
            print "Generate Data Set from " + base_list + "\nusing this bed file : " + bedFile
            if args.testsamplename != None:
                testsamplename = args.testsamplename
                createDataSetFromList_TEST(base_list,bedFile,"1")
                createDataSetFromList_TEST(base_list,bedFile,"2")
                classifying_test()
            else:
                createDataSetFromList(base_list,bedFile)
                classifying()


#    outFileName = args.class_file
#    if args.method == "clustering":
#        print "Classifying data set based on kNN ",str(args.KNN)
#        clustering(int(args.KNN))
#    elif args.method =="classifying":
#        classifying()

 #       if args.PDF_flag != None:
#    output_filter()
    pdf_tag = out_tag
    generate_R_scripts()
    run_R_scripts()
#	remove_internal_files()
