#!/usr/bin/env python

### This is a program to lay out the lattice ###
### geometry for the system                  ###

import numpy as np

def lattice_map(L,nbr):

    SZE = L*L

    NBR_MAX = len(nbr[0])

    #Automated assignment of 6 nbrs and 6
    #nearest nbrs for each lattice point
    for i in range(0,SZE):

        if (np.floor((1.0*i-1.0)/L)==np.floor((1.0*i)/L)):
            nbr[i,0]=i-1
        elif (np.floor((1.0*i-1.0)/L)!=np.floor((1.0*i)/L)):
            nbr[i,0]=i-1+L

        if (np.floor((1.0*i+1.0)/L)==np.floor((1.0*i)/L)):
            nbr[i,3]=i+1
        elif (np.floor((1.0*i+1.0)/L)!=np.floor((1.0*i)/L)):
            nbr[i,3]=i+1-L


        if ((i+L)>=SZE):
            NUM=i+L-SZE		
            nbr[i,2]=NUM
            if (np.floor((1.0*NUM-1.0)/L)==np.floor((1.0*NUM)/L)):
                nbr[i,1]=NUM-1
            elif (np.floor((1.0*NUM-1.0)/L)!=np.floor((1.0*NUM)/L)):
                nbr[i,1]=NUM-1+L

            if NBR_MAX == 12:

                if (np.floor((1.0*NUM+1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+2]=NUM+1
                elif (np.floor((1.0*NUM+1.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+2]=NUM+1-L
                if (np.floor((1.0*NUM-2.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+0]=NUM-2
                elif (np.floor((1.0*NUM-2.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+0]=NUM-2+L

        elif ((i+L)<SZE):
            NUM=i+L		
            nbr[i,2]=NUM
            if (np.floor((1.0*NUM-1.0)/L)==np.floor((1.0*NUM)/L)):
                nbr[i,1]=NUM-1
            elif (np.floor((1.0*NUM-1.0)/L)!=np.floor((1.0*NUM)/L)):
                nbr[i,1]=NUM-1+L
            
            if NBR_MAX == 12:

                if (np.floor((1.0*NUM+1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+2]=NUM+1
                elif (np.floor((1.0*NUM+1.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+2]=NUM+1-L
                if (np.floor((1.0*NUM-2.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+0]=NUM-2
                elif (np.floor((1.0*NUM-2.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+0]=NUM-2+L

        if ((i-L)<0):
            NUM=i-L+SZE		
            nbr[i,5]=NUM
            if (np.floor((1.0*NUM+1.0)/L)==np.floor((1.0*NUM)/L)):
                nbr[i,4]=NUM+1
            elif (np.floor((1.0*NUM+1.0)/L)!=np.floor((1.0*NUM)/L)):
                nbr[i,4]=NUM+1-L

            if NBR_MAX == 12:

                if (np.floor((1.0*NUM+2.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+3]=NUM+2
                elif (np.floor((1.0*NUM+2.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+3]=NUM+2-L
                if (np.floor((1.0*NUM-1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+5]=NUM-1
                elif (np.floor((1.0*NUM-1.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+5]=NUM-1+L

        elif ((i-L)>=0):
            NUM=i-L		
            nbr[i,5]=NUM
            if (np.floor((1.0*NUM+1.0)/L)==np.floor((1.0*NUM)/L)):
                nbr[i,4]=NUM+1
            elif (np.floor((1.0*NUM+1.0)/L)!=np.floor((1.0*NUM)/L)):
                nbr[i,4]=NUM+1-L

            if NBR_MAX == 12:

                if (np.floor((1.0*NUM+2.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+3]=NUM+2
                elif (np.floor((1.0*NUM+2.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+3]=NUM+2-L
                if (np.floor((1.0*NUM-1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+5]=NUM-1
                elif (np.floor((1.0*NUM-1.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+5]=NUM-1+L

        if NBR_MAX == 12:

            if ((i+2*L)>=SZE):
                NUM=i+2*L-SZE		
                if (np.floor((1.0*NUM-1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+1]=NUM-1
                elif (np.floor((1.0*NUM-3.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+1]=NUM-1+L

            elif ((i+2*L)<SZE):
                NUM=i+2*L		
                if (np.floor((1.0*NUM-1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+1]=NUM-1
                elif (np.floor((1.0*NUM-3.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+1]=NUM-1+L

            if ((i-2*L)<0):
                NUM=i-2*L+L*L		
                if (np.floor((1.0*NUM+1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+4]=NUM+1
                elif (np.floor((1.0*NUM+1.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+4]=NUM+1-L

            elif ((i-2*L)>=0):
                NUM=i-2*L		
                if (np.floor((1.0*NUM+1.0)/L)==np.floor((1.0*NUM)/L)):
                    nbr[i,6+4]=NUM+1
                elif (np.floor((1.0*NUM+1.0)/L)!=np.floor((1.0*NUM)/L)):
                    nbr[i,6+4]=NUM+1-L

