#cython: boundscheck=False, wraparound=False, nonecheck=False
#cython: cdivision=True, profile=True

### This is a program to lay out the lattice ###
### geometry for the system                  ###

from libc.math cimport floor

def design(long LSYS,long[:,:] nbr):

    #declaring C datatypes for used variables
    cdef long i, NSYS, NUM
    NSYS = LSYS*LSYS

    #Automated assignment of 6 nbrs and 6
    #nearest nbrs for each lattice point
    for i in range(0,NSYS):

        if (floor((1.0*i-1.0)/LSYS)==floor((1.0*i)/LSYS)):
            nbr[i,0]=i-1
        elif (floor((1.0*i-1.0)/LSYS)!=floor((1.0*i)/LSYS)):
            nbr[i,0]=i-1+LSYS

        if (floor((1.0*i+1.0)/LSYS)==floor((1.0*i)/LSYS)):
            nbr[i,1]=i+1
        elif (floor((1.0*i+1.0)/LSYS)!=floor((1.0*i)/LSYS)):
            nbr[i,1]=i+1-LSYS


        if ((i+LSYS)>=NSYS):
            NUM=i+LSYS-LSYS*LSYS		
            nbr[i,2]=NUM
            if (floor((1.0*NUM-1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,3]=NUM-1
            elif (floor((1.0*NUM-1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,3]=NUM-1+LSYS

            if (floor((1.0*NUM+1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,6+0]=NUM+1
            elif (floor((1.0*NUM+1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,6+0]=NUM+1-LSYS
            if (floor((1.0*NUM-2.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,7]=NUM-2
            elif (floor((1.0*NUM-2.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,7]=NUM-2+LSYS

        elif ((i+LSYS)<NSYS):
            NUM=i+LSYS		
            nbr[i,2]=NUM
            if (floor((1.0*NUM-1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,3]=NUM-1
            elif (floor((1.0*NUM-1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,3]=NUM-1+LSYS

            if (floor((1.0*NUM+1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,6+0]=NUM+1
            elif (floor((1.0*NUM+1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,6+0]=NUM+1-LSYS
            if (floor((1.0*NUM-2.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,7]=NUM-2
            elif (floor((1.0*NUM-2.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,7]=NUM-2+LSYS

        if ((i-LSYS)<0):
            NUM=i-LSYS+LSYS*LSYS		
            nbr[i,4]=NUM
            if (floor((1.0*NUM+1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,5]=NUM+1
            elif (floor((1.0*NUM+1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,5]=NUM+1-LSYS

            if (floor((1.0*NUM+2.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,8]=NUM+2
            elif (floor((1.0*NUM+2.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,8]=NUM+2-LSYS
            if (floor((1.0*NUM-1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,9]=NUM-1
            elif (floor((1.0*NUM-1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,9]=NUM-1+LSYS

        elif ((i-LSYS)>=0):
            NUM=i-LSYS		
            nbr[i,4]=NUM
            if (floor((1.0*NUM+1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,5]=NUM+1
            elif (floor((1.0*NUM+1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,5]=NUM+1-LSYS

            if (floor((1.0*NUM+2.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,8]=NUM+2
            elif (floor((1.0*NUM+2.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,8]=NUM+2-LSYS
            if (floor((1.0*NUM-1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,9]=NUM-1
            elif (floor((1.0*NUM-1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,9]=NUM-1+LSYS

        if ((i+2*LSYS)>=NSYS):
            NUM=i+2*LSYS-LSYS*LSYS		
            if (floor((1.0*NUM-1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,10]=NUM-1
            elif (floor((1.0*NUM-3.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,10]=NUM-1+LSYS

        elif ((i+2*LSYS)<NSYS):
            NUM=i+2*LSYS		
            if (floor((1.0*NUM-1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,10]=NUM-1
            elif (floor((1.0*NUM-3.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,10]=NUM-1+LSYS

        if ((i-2*LSYS)<0):
            NUM=i-2*LSYS+LSYS*LSYS		
            if (floor((1.0*NUM+1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,11]=NUM+1
            elif (floor((1.0*NUM+1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,11]=NUM+1-LSYS

        elif ((i-2*LSYS)>=0):
            NUM=i-2*LSYS		
            if (floor((1.0*NUM+1.0)/LSYS)==floor((1.0*NUM)/LSYS)):
                nbr[i,11]=NUM+1
            elif (floor((1.0*NUM+1.0)/LSYS)!=floor((1.0*NUM)/LSYS)):
                nbr[i,11]=NUM+1-LSYS


