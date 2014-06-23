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

    def get_stat(self, mpireduce=None, tagofs=0):
        if (self.sw==0):
            return
        if (mpireduce):
            X=self.x*1.0
            XX=self.xx*1.0
            rank=mpireduce.Get_rank()
            size=mpireduce.Get_size()
            if rank == 0:
                for i in range(1,size):
                    data=X*0.0
                    mpireduce.Recv(data, source=i, tag=tagofs+1)
                    X+=data
                    data=XX*0.0
                    mpireduce.Recv(data, source=i, tag=tagofs+2)
                    XX+=data
            else:
                mpireduce.Send(X,dest=0,tag=tagofs+1)
                mpireduce.Send(XX,dest=0,tag=tagofs+2)
                ## need to get mean, error back
                data=X*0.0
                mpireduce.Recv(data,source=0, tag=tagofs+3)
                self.mean=data
                data=X*0.0
                mpireduce.Recv(data,source=0, tag=tagofs+4)
                self.err=data

                return

            
            SW=self.sw*size
        else:
            X=self.x*1.0
            XX=self.xx*1.0
            SW=self.sw

        m=X/SW
        mm=XX/SW
        mm-=m*m
        self.mean=m
        self.err=sqrt(mm/SW)
        if (mpireduce):
            for i in range(1,size):
                mpireduce.Send(self.mean,dest=i,tag=tagofs+3)
                mpireduce.Send(self.err,dest=i,tag=tagofs+4)

                


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

        
    
 
