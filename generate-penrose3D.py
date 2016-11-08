import math
import csv
import operator
import time
from scipy.optimize import linprog

# One parameter:
len = 5    # Size of the bounding box

output_filename = "penrose3d-5x5x5.csv"


r = (math.sqrt(5)+1)/2
q = (math.sqrt(5)-1)/2

def projection( p ):
    a = math.sqrt(2/(5+math.sqrt(5)))
    return [a*(p[0]+r*p[1]-p[3]+r*p[4]),\
            a*(r*p[0]+p[2]+r*p[3]-p[5]),\
            a*(p[1]+r*p[2]-p[4]+r*p[5])]    

def proj_in_box( p ):
    tp = projection(p)
    return tp[0] >= 0 and tp[1] >= 0 and tp[2] >= 0 and\
            tp[0] <= len and tp[1] <= len and tp[2] <= len

def add_vec(vec1, vec2):
    return tuple(map(operator.add, vec1, vec2))

def Generate3DPenrose():
    # Next 5 parameters (reals) define the 3D space into which penrose is projected.
    # Varying them will change the geometry of the resulting tiling.
    g0 = 0.1
    g1 = 0.2
    g2 = 0.5
    g3 = 0.15
    g4 = 0.3
    
    csvfile = open('output_filename', 'wb')
    csvwriter = csv.writer(csvfile)
    p0min = int(math.floor(0.276*g0-0.447*(g1+g2+g3+g4)))
    p1min = int(math.floor(0.724*g1+0.447*g3))
    p2min = int(math.floor(0.724*g2+0.447*g4))
    p0max = int(math.ceil(0.688*len+0.276*g0-0.447*(g1+g2+g3+g4)))
    p1max = int(math.ceil(0.688*len+0.724*g1+0.447*g3))
    p2max = int(math.ceil(0.688*len+0.724*g2+0.447*g4))
    
    vertices = dict()
    for p0 in range(p0min,p0max):
        for p1 in range(p1min,p1max):
            for p2 in range(p2min,p2max):
                p3_ = math.floor(q*p0-p1+q*p2-q*g0+g1-q*g2+g3)
                p4_ = math.floor(q*p0+q*p1-p2-q*g0-q*g1+g2+g4)
                p5_ = math.floor(-p0+q*p1+q*p2-r*g1-r*g2-g3-g4)
                for i in range(-3,4):
                    for j in range(-3,4):
                        for k in range(-3,4):
                            p3 = p3_ + i
                            p4 = p4_ + j
                            p5 = p5_ + k
                            if proj_in_box([p0,p1,p2,p3,p4,p5]):
                                c = [0, 0, 0]
                                A = [[q, -1, q], [q, q, -1], [-1, q, q], [-q, 1, -q], [-q, -q, 1], [1, -q, -q]]
                                b = [p3+1-(-q*g0+g1-q*g2+g3), p4+1-(-q*g0-q*g1+g2+g4), p5+1-(-r*g1-r*g2-g3-g4),\
                                     -p3-q*g0+g1-q*g2+g3, -p4-q*g0-q*g1+g2+g4, -p5-r*g1-r*g2-g3-g4]
                                
                                x0_bounds = (p0, p0+1)
                                x1_bounds = (p1, p1+1)
                                x2_bounds = (p2, p2+1)
                                
                                res = linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds, x2_bounds))
                                if res.success:
                                    vertices[(p0, p1, p2, p3, p4, p5)] = projection((p0, p1, p2, p3, p4, p5))

    edges = ((1,0,0,0,0,0), (0,1,0,0,0,0), (0,0,1,0,0,0), (0,0,0,1,0,0), (0,0,0,0,1,0), (0,0,0,0,0,1))
    for v, pr0 in vertices.iteritems():
        for e in edges:
            vi = add_vec(v, e)
            if vi in vertices:
                pri = vertices[vi]
                csvwriter.writerow([pr0[0], pr0[1], pr0[2], pri[0], pri[1], pri[2]])

start_time = time.time()
Generate3DPenrose()
print("--- %s seconds ---" % (time.time() - start_time))
