import os
import random
from math import exp, pow, sqrt, fabs
import numpy as np
import ShaKer.simushape as ss
from ShaKer.rna_tools.rnasubopt import  rnasubopt

# hacked version of sukosds SHAPE simulation method
# http://users-birc.au.dk/zs/SHAPEsimulations/ 

def expCDF(x):
	#Exponential distribution CDF
	lamb = 0.681211; #lambda
	dist = 1-exp(-x/lamb);
	return dist; #minus desired value so we can seek minimum

def gevCDFouter(x):
	#Generalized Extreme Value distribution CDF, outer pairs
	xi = 0.821235;
	oneoverxi = 1/xi;
	sigma = 0.113916;
	mu = 0.0901397;
	dist = exp(-1*pow(1 + xi*(x-mu)/sigma,-oneoverxi) ); 
	return dist; #minus desired value so we can seek minimum

def gevCDFinner(x):
	#Generalized Extreme Value distribution CDF, inner pairs
	xi = 0.762581;
	oneoverxi = 1/xi;
	sigma = 0.0492536;
	mu = 0.0395857;
	dist = exp(-1*pow(1 + xi*(x-mu)/sigma,-oneoverxi) ); 
	return dist; #minus desired value so we can seek minimum
	
phi = (1 + sqrt(5)) / 2;
resphi = 2 - phi;
 
# x1 and x3 are the current bounds; the minimum is between them.
# x2 is the center point, which is closer to x1 than to x3
def goldenSectionSearch(f, desired, x1, x2, x3, tau):
	#calculate new potential center point
	x4 = x2 + resphi * (x3 - x2);
	if fabs(x3 - x1) < tau * (fabs(x2) + fabs(x4)):
		return (x3 + x1) / 2;
	if fabs(f(x4)-desired) < fabs(f(x2)-desired):
		return goldenSectionSearch(f, desired, x2, x4, x3, tau);
	else:
		return goldenSectionSearch(f, desired, x4, x2, x1, tau)

def randomSHAPEinnerpairing():
	random.seed();
	randomnumber = random.random();	
	return goldenSectionSearch(gevCDFinner, randomnumber, 0, resphi*10, 10, sqrt(1e-10)); 
	#return 0.1;

def randomSHAPEouterpairing():
	random.seed();
	randomnumber = random.random();	
	return goldenSectionSearch(gevCDFouter, randomnumber, 0, resphi*10, 10, sqrt(1e-10)); 
	#return 0.1;
	
def randomSHAPEunpaired():
	random.seed();
	randomnumber = random.random();	
	return goldenSectionSearch(expCDF, randomnumber, 0, resphi*10, 10, sqrt(1e-10)); 
	#return 10;
	
def generateValue(n,m,o):
	if(n>0):
		if(m>0 and o>0):
		#inner pairing
			return randomSHAPEinnerpairing();
		else:
		#outer pairing
			return randomSHAPEouterpairing();
	else:
		#not pairing
		return randomSHAPEunpaired();




def db_to_suko(db):
	#encoding looks like this [0, 20, 19, 18, 17, 16, 15, 14, 13, 0, 0, 0, 9, 8, 7, 6, 5, 4, 3, 2]
        result = []

        stack = []
        stack2=[]
        for i,e in enumerate(db):
            if e == '.':
                result.append(0)
            if e == '(':
                result.append(-1)
                stack.append(i)
            if e == '[':
                result.append(-1)
                stack2.append(i)
            if e == ')':
                otherid = stack.pop()
                result[otherid] = i+1
                result.append(otherid+1)
            if e == ']':
                otherid = stack2.pop()
                result[otherid] = i + 1
                result.append(otherid + 1)
        return result
	


def sukosd(dotbracket):
    real_pairing  = db_to_suko(dotbracket)
    length = len(real_pairing);
    data = []
    for i in range(0,length):
            if i > 0:
                    if real_pairing[i-1] == real_pairing[i]+1:
                            m = 1;
                    else:
                            m = 0;
            else:
                    m = 0;
            if i < length-1:
                    if real_pairing[i+1] == real_pairing[i]-1:
                            o = 1;
                    else:
                            o = 0;
            else:
                    o = 0;
            data.append( generateValue(real_pairing[i],m,o))
    return np.array(data)


def predict_Suko(sequence,seq_to_db_function= rnasubopt):
    db_list = seq_to_db_function(sequence)
    struct_proba = ss.probabilities_of_structures(sequence, db_list)
    structures, weights = zip(*struct_proba)
    shapes_Suko = [sukosd(x) for x in db_list]
    return ss.weighted_average(weights, shapes_Suko)



def predict_tmp(sequence,seq_to_db_function= rnasubopt):
    db_list = seq_to_db_function(sequence)
    struct_proba = ss.probabilities_of_structures(sequence, db_list)
    structures, weights = zip(*struct_proba)
    shapes = [ np.array([ 1.0 if chr == "." else 0.0 for chr in x ]) for x in db_list]
    return ss.weighted_average(weights, shapes)



