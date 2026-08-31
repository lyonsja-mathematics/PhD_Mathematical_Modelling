#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def trapezoidal_rule_integration(array, step):
    
    import numpy as np
    if np.sqrt(len(array))%1 != 0:
    
        raise Exception("Array must be a square matrix.")
        
    else:
    
        Nx = int(np.sqrt(len(array)))
        array=array.reshape(Nx,Nx).T
    
        return step*step/4*(array[0,0] + array[Nx-1,0] + array[0,Nx-1] + array[Nx-1,Nx-1] + np.sum(2*array[0,1:(Nx-1)]) + np.sum(2*array[1:(Nx-1),0]) + np.sum(2*array[Nx-1,1:(Nx-1)]) + np.sum(2*array[1:(Nx-1),Nx-1]) + np.sum(4*array[1:(Nx-1),1:(Nx-1)]))

def simpsons_rule_1D(array_1D, step):
    
    intpts=list(range(1,len(array_1D)-1))
    times_4=[i for i in intpts if i%2!=0]
    times_2=[i for i in intpts if i%2==0]
    
    return step/3*(array_1D[0] + 4*sum(array_1D[times_4]) + 2*sum(array_1D[times_2]) + array_1D[len(array_1D)-1])
