import math
import csv
import operator
import time
from scipy.optimize import linprog

# One parameter:
len = 10    # Size of the bounding box

output_filename = "penrose2d.csv"

def projection( p ):
   return [p[0]+p[1]*math.cos(2*math.pi/5)+p[2]*math.cos(4*math.pi/5)+p[3]*math.cos(6*math.pi/5)+p[4]*math.cos(8*math.pi/5),\
           p[1]*math.sin(2*math.pi/5)+p[2]*math.sin(4*math.pi/5)+p[3]*math.sin(6*math.pi/5)+p[4]*math.sin(8*math.pi/5)]

def proj_in_box( p, len ):
    tp = projection(p)
    return tp[0] >= 0 and tp[1] >= 0 and\
        tp[0] <= len and tp[1] <= len

def add_vec(vec1, vec2):
    return tuple(map(operator.add, vec1, vec2))

def Generate2DPenrose():
    r = (math.sqrt(5)+1)/2
    q = (math.sqrt(5)-1)/2

    # Next 4 parameters (reals) define the plane into which penrose is projected.
    # Varying them will change the geometry of the resulting tiling.
    g0 = 0.1
    g1 = 0.2
    g2 = 0.5
    g3 = 0.15

    csvfile = open(output_filename, 'wb')
    csvwriter = csv.writer(csvfile)
    
    p0min = int(math.floor(0.724*g0+0.447*g2+0.447*g3))
    p0max = int(math.ceil(0.4*len+0.724*g0+0.447*g2+0.447*g3))
    p1min = int(math.floor(-0.447*g0+0.276*g1-0.447*g2-0.447*10**(-16)*g3))
    p1max = int(math.ceil(0.504*len-0.447*g0+0.276*g1-0.447*g2-0.447*10**(-16)*g3))
    
    vertices = dict()
    for p0 in range(p0min,p0max):
        for p1 in range(p1min,p1max):
            p2_ = math.floor(-p0+q*p1+g0-q*g1+g2)
            p3_ = math.floor(-q*p0-q*p1+q*g0+q*g1+g3)
            p4_ = math.floor(q*p0-p1-r*g0-g2-g3)
            for i in range(-3,4):
                for j in range(-3,4):
                    for k in range(-3,4):
                        p2 = p2_ + i
                        p3 = p3_ + j
                        p4 = p4_ + k
                        if proj_in_box([p0,p1,p2,p3,p4], len):
                            c = [1, 1]
                            A = [[-1, q], [-q, -q], [q, -1],\
                                 [1, -q], [q, q], [-q, 1]]
                            b = [p2+1-g0+q*g1-g2, p3+1-q*g0-q*g1-g3, p4+1+r*g0+g2+g3,\
                                 -(p2-g0+q*g1-g2), -(p3-q*g0-q*g1-g3), -(p4+r*g0+g2+g3)]

                            x0_bounds = (p0, p0+1)
                            x1_bounds = (p1, p1+1)

                            res = linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds))
                            if res.success:
                                vertices[(p0, p1, p2, p3, p4)] = projection((p0, p1, p2, p3, p4))

    edges = ((1,0,0,0,0), (0,1,0,0,0), (0,0,1,0,0), (0,0,0,1,0), (0,0,0,0,1))
    for v, pr0 in vertices.iteritems():
        for e in edges:
            vi = add_vec(v, e)
            if vi in vertices:
                pri = vertices[vi]
                csvwriter.writerow([pr0[0], pr0[1], pri[0], pri[1]])

start_time = time.time()
Generate2DPenrose()
print("--- %s seconds ---" % (time.time() - start_time))
