#!/usr/bin/python3
# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20210603
Description:
This Script of nanopore data analysis pipeline is developed by Legendzdy.
"""
#%%
import argparse, time, os, re, sys
import subprocess

class unifiedInterface(object):
    def __init__(self, anySoft, fowardDir):
        self.fowardDir = fowardDir
        self.anySoft = anySoft
    
    def getWorkDir(self):
        with open('/inputs/config.csv', 'r') as file:
            lines = file.readlines()
            header = lines[0].strip().split(',') 
            prefixIndex = header.index('prefix')
            WorkDir = "OUT_" + lines[1].strip().split(',')[prefixIndex]
        return WorkDir
    
    def check(self):
        workDir=self.getWorkDir()
        workDirPath = "/outputs" + "/" + workDir
        os.chdir(workDirPath)
        dirList = os.listdir(os.getcwd())
        if re.search(r"results_", self.fowardDir) != None:
            fowardDirString= self.fowardDir
            fowardDirRe= re.compile(fowardDirString)
            fowardPath = os.path.join(os.getcwd(), fowardDirString)
        else:
            fowardDirString= "_" + self.fowardDir + "_"
            fowardDirRe= re.compile(fowardDirString)
        
        if self.fowardDir == "input":
            outNum = "01"
            outDir = "results_" + outNum + "_" + self.anySoft
            if not any(re.search(outDir, i) is not None for i in dirList):
                os.mkdir(outDir)
                fowardPath = os.path.join(os.getcwd(), "input")
                targetPath = os.path.join(os.getcwd(), outDir)
                os.chdir(targetPath)
                print("OUTPUT: {}".format(outDir))
            else:
                fowardPath = "exist"
        else:
            if any(fowardDirRe.search(i) is not None for i in dirList):
                for i in dirList:
                    if fowardDirRe.search(i) != None:
                        if re.search(r"_success", i) != None:
                            nowNum = re.findall(r"_(\d+)_", i)[0]
                            outNum = self.twoDigits(int(nowNum) + 1)
                            outDir = "results_" + outNum + "_" + self.anySoft
                            if not any(re.search(outDir, i) is not None for i in dirList):
                                os.mkdir(outDir)
                                fowardDirNew = i
                                fowardPath = os.path.join(os.getcwd(), fowardDirNew)
                                targetPath = os.path.join(os.getcwd(), outDir)
                                os.chdir(targetPath)
                                print("OUTPUT: {}".format(outDir))
                            else:
                                fowardPath = "exist"
                        else:
                            print("The last step is failed!")
                            sys.exit(1)
            else:
                raise Exception("The input folder is not exist!")
        return fowardPath

    def twoDigits(self, num):
        if num < 10:
            return "0" + str(num)
        else:
            return str(num)
    
    def run(self):
        start = time.time()
        fowardPath = self.check()
        targetPath = os.getcwd()
        if fowardPath == "exist":
            targetDir=os.path.basename(targetPath)
            print("{} is exist!".format(targetDir))
        else:
            try:
                command = ['bash','/root/legendzdy_module.sh', '-i', fowardPath, '-o', targetPath]
                subprocess.run(command, check=True)
                end = time.time()
                print("Total time: {:.2f} min".format((end - start) / 60))
                timeFile = self.anySoft + "_run_time.tsv"
                with open(timeFile, 'w') as f:
                    f.write("Total time: {:.2f} min".format((end - start) / 60))
                targetDir=os.getcwd()
                os.chdir("..")
                os.rename(targetDir, targetDir + "_success")
            except subprocess.CalledProcessError as e:
                print(f"An error occurred: {e}")
                print(f"Error output: {e.stderr}")
                targetDir=os.getcwd()
                os.chdir("..")
                os.rename(targetDir, targetDir + "_failed")
                sys.exit(1)

class softDemo(unifiedInterface):
    def __init__(self, anySoft, fowardDir):
        super().__init__(anySoft, fowardDir)

    def demo(self):
        super().run()
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='clean raw data with cutadapt.'
    )
    parser.add_argument('--version', action='version', version='cutadapt Version: %(prog)s 4.5')
    parser.add_argument('-f', '--fowardDir',  required=True, default="input", help='The foward dir.')
    parser.add_argument('--anySoft', default="cutadapt", help=argparse.SUPPRESS)
    args = parser.parse_args()
    
    softDemo(args.anySoft, args.fowardDir).demo()