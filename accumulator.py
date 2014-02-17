from scipy import *


class Accumulator:
    def __init__ (self,name,Nch=30):
        self.name=name
        self.sw=0

    def accumulate(self, obj):
        if (self.sw>0):
            self.x+=obj
            self.xx+=obj*obj
        else:
            self.x=obj
            self.xx=obj*obj
        self.sw+=1 

    def get_stat(self):
        if (self.sw==0):
            return
        self.x/=self.sw
        self.xx/=self.sw
        self.xx-=self.x*self.x
        self.mean=self.x
        self.err=sqrt(self.xx/self.sw)


    def print_stat(self):
        if (self.sw==0):
            return
        self.get_stat()
        print '--------------------------'
        print "Stats for ",self.name
        
        print "Mean:"
        print self.mean
        print "Error:"
        print self.err






class ChunkAccumulator:
    def __init__ (self,name,nch):
        self.name=name
        self.Nch=nch
        self.sw=zeros((nch,))
        self.vals = [None]*nch
        self.cc=0

    def accumulate(self, obj):
        if self.vals[self.cc]==None:
            self.vals[self.cc]=obj
        else:
            self.vals[self.cc]+=obj
        self.sw[self.cc]+=1
        self.cc=(self.cc+1)%self.Nch

    def print_stat(self):
        print '--------------------------'
        print "Stats for ",self.name
        ## first get the total mean
        tot=reduce(lambda x,y:x+y, self.vals)
        tot/=self.sw.sum()
        print "Mean:"
        print tot
        ## now the errors
        vals=[val/we for val,we in zip(self.vals,self.sw)]
        vals-=tot
        valssq=vals*vals
        stddev=reduce(lambda x,y:x+y,valssq)/self.Nch
        err=sqrt(stddev/(self.Nch-1))
        print "Error:"
        print err

        
    
 
