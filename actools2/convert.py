# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 20:24:11 2014

@author: Maxim
"""

def m_to_ft(length):
    return 3.2808399*float(length)

def ft_to_m(length):
    return float(length)/3.2808399

def sqm_to_sqft(area):
    return 10.76391044943201*float(area)

def sqft_to_sqm(area):
    return float(area)/10.76391044943201

def kg_to_lb(mass):
    return 2.20462262 * float(mass)

def lb_to_kg(mass):
    return float(mass)/2.20462262


def run_test1():
    print m_to_ft(1)
    print sqm_to_sqft(15)

if __name__=="__main__":
    run_test1()