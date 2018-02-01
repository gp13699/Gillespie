import random as r
import math as m
import numpy as np
import matplotlib.pyplot as plt
import os

#Make more elegant
base_path = "/Users/gregoryionovichpage/Desktop/" #change to your own path
filenameJ = "ts5p1J.txt"  #absolute path to mathematica outputs for ts=5.1 !!
filenameH = "ts5p1H.txt" 
filenameMsq = "Var5p1.txt"
path_to_fileJ = os.path.join(base_path, filenameJ)
path_to_fileH = os.path.join(base_path, filenameH) 
path_to_fileMsq = os.path.join(base_path, filenameMsq) 

#fJ = open(path_to_fileJ , 'r')#opening file, read mode
#fH = open(path_to_fileH , 'r')#opening file, read mode
#file = open(path_to_file , 'r')

print (path_to_fileJ)
print (path_to_fileH)
print (path_to_fileMsq)

Decide= 'J'
ts=5.1 #stopping time WARNING, if changing, you might also want new mathematica input files to read in.
tau=1 #transit time
S=10000 #number of samples
Lmbd=10 #Poisson intensity upper limit

class ClassHJ(object):

    def __init__(self,Lmbd,tau,ts,S,Decide):
        
        self.Lmbd=Lmbd
        self.tau=tau
        self.ts=ts
        self.S=S
    
    if Decide == 'H':  

        def Channel(self,lmbd):
            #want output, i to be the number of exiting particles GIVEN a blockage occurs before ts
            a=0 #for checking contributions
            i=0#reset to cero
            Ttot=0
            while True: 

                E=-(m.log(r.random()))/lmbd #E for event

                TtotO=Ttot #store old total time    
                Ttot+=E #add to total time the sum of event times
                Diff= Ttot-TtotO #time difference between events

                i+=1 #is number of particles that have entered

                if(Ttot >= self.ts): 
                    #here, need to be careful
                    #saying, if a particle enters at ts or beyond, we have had no blockage
                    #means that condition for h(0,ts) is not satisfied. 
                    #set i to an impossible value.
                    i= -1
                    #returning this value will force simu to go back one sample, and try again.
                    #a=i

                    break

                if(Diff <= self.tau and i !=1): #'normal blocking condition'
                    #it takes two particles to make a blockage, so
                    i-=2 #this number exited

                    break

            Ttot=0
            return i;
        
    if Decide == 'J': 

        def Channel(self,lmbd):
            #want output, i to be the number of exiting particles
            a=0 #for checking contributions
            i=0
            Ttot=0
            while True: 

                E=-(m.log(r.random()))/lmbd #E for event

                TtotO=Ttot #store old total time    
                Ttot+=E #add to total time the sum of event times
                Diff= Ttot-TtotO #time difference between events

                i+=1 #is number of particles that have entered

                if(Diff <= self.tau and i !=1): #'stronger condition'
                    #it takes two particles to make a blockage, so
                    i-=2 #this number exited
                    #a=i

                    break

                if(Ttot >= self.ts-self.tau):
                    #last event needs to be subtracted
                    i-=1
                    #a=i

                    break

            #need a distinciton between number of particles that exited at ts against blockage?

            Ttot=0
            return i;
        
    if Decide == 'G': 
        #creating simu for K particles exiting AND that it has remained open. 
        
        def Channel(self,lmbd):
            #
            a=0 #for checking contributions
            i=0
            Ttot=0
            while True: 

                E=-(m.log(r.random()))/lmbd #E for event

                TtotO=Ttot #store old total time    
                Ttot+=E #add to total time the sum of event times
                Diff= Ttot-TtotO #time difference between events

                i+=1 #is number of particles that have entered
                
                if(Diff <= self.tau and i !=1): #'stronger condition'
                    #if there is a condition for a blockage, reject the sample!
                    i= -1 #set to impossible value

                    break
                    
                    
                if(Ttot >= self.ts-self.tau): #there could still be a blockage in final time interval. ACCOUNT
                    #last event needs to be subtracted as it won't exit.
                    Ttot+=-(m.log(r.random()))/lmbd  #generate last event
                    #a=i
                    if (Ttot >= self.ts):
                        i-=1
                    else:
                        i= -1
                    #need to generate one more event. if last event is before ts, send -1.
                    #is last event is after ts, do i-=1 
                    
                    break
                    


            #need a distinciton between number of particles that exited at ts against blockage?

            Ttot=0
            return i;
        
    
    def Sample(self,lmbd,PSum,PSumSq):#need new arguments for check

        for j in range(self.S-1):

            P = self.Channel(lmbd)#no. of particles avant block, OR impossible value.

            #print(P,j)
            
            if(P == -1):#if impossible value (define as qqch autre?)
                j-=1 #forces the sanple to be taken again.

            else:#continue as normal, include new values to sums

                PSum += P
                PSumSq += P*P 
                             
        PAv=PSum/self.S
        PAvSq=PSumSq/self.S
        Var= PAvSq - (PAv*PAv) #variance of <m>

        return PAvSq, Var
    
    def Iteration(self):

        #variables to modify are lmbd,
        lmbd=.0001 #lambda iterator 
        K=200 #for splitting up data for representation later
        
        x=[]#empty lists for plotting later
        y=[]
        z=[]
        X=[]

        dataMsq = np.genfromtxt(path_to_fileMsq)#Read from txt mathematica output
        #dataJ = np.genfromtxt(path_to_fileJ)#Read from txt mathematica output
        xM=dataMsq[:,0]
        yM=dataMsq[:,1] 
        #xM1=dataJ[:,0]
        #yM1=dataJ[:,1]
        for k in range(K):

            PSum=0
            PSumSq=0
            lmbd+=(self.Lmbd/K)

            self.Channel(lmbd) 
            A,B = self.Sample(lmbd,PSum, PSumSq)

            x.append(lmbd)
            X.append( (1+np.exp(lmbd))/( (1-np.exp(lmbd)) * (1-np.exp(lmbd)) ) ) # limit output for high t
            
            y.append(A)
            z.append(B)
            



            #f.write(str(lmbd))#Output text file, finish. 
            
        plt.subplot(221) #look up meaning of number
        plt.grid(True)
        
        
        plt.plot(x, z)
        #plt.plot(x, yM,'r:')
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.title('log')        
        plt.plot(xM, yM,'r:')#plot analytical outputs
        #plt.plot(xM1, yM1,'g:')#plot analytical outputs
        #plt.errorbar(x, y,np.sqrt(z))
        plt.ylabel('<Var>')
        plt.xlabel('Lambda')
        plt.show()
        #fH.close()


ClassHJ(Lmbd,tau,ts,S,Decide).Iteration()