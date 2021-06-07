from mpi4py import MPI
import numpy
import sys
import json
WORKTAG = 0
DIETAG = 1

class Work():
    '''
    pop out job from the bottom
    '''
    def __init__(self, work_items):
        '''
        work_items must be array- or list-like
        '''
        try:
            self.work_items = work_items.tolist()
        except:
            self.work_items = work_items

    def get_next_item(self):
        if len(self.work_items) == 0:
            return None
        return self.work_items.pop()

def master(wi):
    all_data = []
    size = MPI.COMM_WORLD.Get_size()
    current_work = Work(wi) 
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    nsend,nrece = 0,0

    for i in range(1, size): 
        anext = current_work.get_next_item()
        if not anext:
            break
        comm.send(obj=anext,
                  dest=i, tag=WORKTAG)
        nsend += 1

    while 1:
        anext = current_work.get_next_item()
        if not anext: break
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        print('data_recv:', data)
        nrece += 1
        #data['x'] = dat [{'index': i, 'x': pi, 'costf': func(np.array(pi))} ]
        all_data.append({'index': len(all_data),'x': data[0], 'costf': data[1]})
        comm.send(obj=anext, dest=status.Get_source(), tag=WORKTAG)
        nsend += 1
    #print 'master sent out all data'

    for i in range(1,size):
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
        #print('data_reciev2:', data)
        nrece += 1
        all_data.append({'index': len(all_data),'x': data[0], 'costf': data[1]})

    #print "# of sent: %i"%nsend
    #print "# of rece: %i"%nrece
    #print 'master received all data'

    for i in range(1,size):
        comm.send(obj=None, dest=i, tag=DIETAG)

    return all_data


def slave(f,*args):
    comm = MPI.COMM_WORLD
    status = MPI.Status()
    while 1:
        data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        # print('data_slave_recv:',data)
        if status.Get_tag():
            break
#        comm.send(obj=f(data,*args), dest=0)
        data.append(f(data, *args))
        comm.send(obj=data, dest=0)

def run(wi,f,*args):
    '''
    wi:   work_items
    f:    the function needs to be implemented on work_items
    args: argument list for f
    '''
    rank = MPI.COMM_WORLD.Get_rank()
    size = MPI.COMM_WORLD.Get_size()
    if size < 2:
        MPI.Finalize()
        sys.exit(-1)

    all_data = None        
    if rank == 0:
        all_data = master(wi)
        with open('data_optimize_grid_scan.txt', 'w') as f:
            json.dump(all_data, f, indent=2) # save the sample points to data file.
        #numpy.savetxt('data_optimize_grid_scan.txt',all_data)   # save the sample points to data file.
    else:
        slave(f,*args)
    MPI.Finalize()
    return all_data

#import numpy as np
#if __name__ == "__main__":

#    def f(data*args):
#        return np.sum(data)+args[0]
#    data = np.ones((10,4))
#    print('data0:',data)
#    run(data,f,3)
#    print('data1:',data)
