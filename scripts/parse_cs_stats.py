#!/usr/bin/env python3
import argparse

def parse_cs(cs):
    m=sub=ins=dele=0
    i=0
    while i<len(cs):
        c=cs[i]
        if c==':':
            i+=1; j=i
            while j<len(cs) and cs[j].isdigit(): j+=1
            m+=int(cs[i:j]); i=j
        elif c=='*':
            sub+=1; i+=3
        elif c=='+':
            i+=1; j=i
            while j<len(cs) and cs[j] not in ":+-*": j+=1
            ins+=(j-i); i=j
        elif c=='-':
            i+=1; j=i
            while j<len(cs) and cs[j] not in ":+-*": j+=1
            dele+=(j-i); i=j
        else:
            i+=1
    return m,sub,ins,dele

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("-o","--out", required=True)
    args=ap.parse_args()

    sub=ins=dele=match=0
    with open(args.paf) as f:
        for line in f:
            for x in line.split("\t"):
                if x.startswith("cs:Z:"):
                    m,s,i,d=parse_cs(x[5:])
                    match+=m; sub+=s; ins+=i; dele+=d

    with open(args.out,"w") as o:
        o.write(f"sub\t{sub}\nins\t{ins}\ndel\t{dele}\n")

if __name__=="__main__":
    main()

