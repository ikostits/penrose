import numpy as np
import math
import time

# Default unit distance between the vertices of the tiling is 1mm
# Use `scale` to increase the length of struts
scale = 6.0

# All the following values will be multiplied by `scale`
t = 0.1                 # Diameter of the inscribed circle in the cross-section 
                        # of the strut (as a fraction of strut length)
                        # = strut width
trim = 2                # Trim the `trim` length off the 
                        # sides of the tiling
plate_thickness = 1     # If value greater than 0, then add two plates on the 
                        # bottom and on the top of the tiling.
                        # Set the value to 0 to get the tiling without the 
                        # bounding plates

input_filename = "penrose3d-15x15x15.csv"
output_filename = "thick-penrose3d-11x11x11-6mm-0.6mm-cap.stl"

precision = 5           # Set the precision to be large enough to distinguish 
                        # between different vertices, but small enough to avoid
                        # the rounding issues. Default value is 5.
eps = 10**(-precision)

class Icosahedron:
    'Vertices are the centers of the dodecahedron-joint'
    __cosA = 1/math.sqrt(5)
    __joint = dict()

    def __init__(self, dirA, dirB, r):
        tau = ((1+math.sqrt(5))/math.sqrt(7)) * r / math.sqrt(dirA.dot(dirA))
        self.vecA = tau*np.asarray(dirA)
        self.vecB = tau*np.asarray(dirB)
        self.vecC = np.dot( rotation_matrix(self.vecA, 2.0*math.pi/5.0), self.vecB )
        self.vecD = np.dot( rotation_matrix(self.vecA, 4.0*math.pi/5.0), self.vecB )
        self.vecE = np.dot( rotation_matrix(self.vecA, 6.0*math.pi/5.0), self.vecB )
        self.vecF = np.dot( rotation_matrix(self.vecA, 8.0*math.pi/5.0), self.vecB )
        self.dodec_verts = [(self.vecA+self.vecB+self.vecC)/(1+2*self.__cosA),\
                            (self.vecA+self.vecC+self.vecD)/(1+2*self.__cosA),\
                            (self.vecA+self.vecD+self.vecE)/(1+2*self.__cosA),\
                            (self.vecA+self.vecE+self.vecF)/(1+2*self.__cosA),\
                            (self.vecA+self.vecF+self.vecB)/(1+2*self.__cosA),\
                            (self.vecB+self.vecC-self.vecE)/(1+2*self.__cosA),\
                            (self.vecC+self.vecD-self.vecF)/(1+2*self.__cosA),\
                            (self.vecD+self.vecE-self.vecB)/(1+2*self.__cosA),\
                            (self.vecE+self.vecF-self.vecC)/(1+2*self.__cosA),\
                            (self.vecF+self.vecB-self.vecD)/(1+2*self.__cosA),\
                            -(self.vecA+self.vecB+self.vecC)/(1+2*self.__cosA),\
                            -(self.vecA+self.vecC+self.vecD)/(1+2*self.__cosA),\
                            -(self.vecA+self.vecD+self.vecE)/(1+2*self.__cosA),\
                            -(self.vecA+self.vecE+self.vecF)/(1+2*self.__cosA),\
                            -(self.vecA+self.vecF+self.vecB)/(1+2*self.__cosA),\
                            -(self.vecB+self.vecC-self.vecE)/(1+2*self.__cosA),\
                            -(self.vecC+self.vecD-self.vecF)/(1+2*self.__cosA),\
                            -(self.vecD+self.vecE-self.vecB)/(1+2*self.__cosA),\
                            -(self.vecE+self.vecF-self.vecC)/(1+2*self.__cosA),\
                            -(self.vecF+self.vecB-self.vecD)/(1+2*self.__cosA)]
        self.map_decagons = dict()
        self._compute_faces()

    def verts(self):
        return [self.vecA, self.vecB, self.vecC, self.vecD, self.vecE, self.vecF,\
                -self.vecA, -self.vecB, -self.vecC, -self.vecD, -self.vecE, -self.vecF]

    def getDirection(self, vec):
        for v in self.verts():
            if equals(vec.dot(v), math.sqrt(vec.dot(vec)*v.dot(v))):
                return v
        return None

    def getDecagon(self, vec):
        key = coords_tuple(vec)
        return self.map_decagons.get(key)

#    def getTriangle(self, vec0, vec1, vec2):
#        key = coords_tuple((vec0+vec1+vec2)/(1+2*self.__cosA))
#        if key in self.map_triangles:
#            return self.map_triangles[key]
#        else:
#            return None

    def _compute_faces(self):
        self.map_decagons[coords_tuple(self.vecA)] = self._getDecagon(self.vecA, self.vecB)
        self.map_decagons[coords_tuple(self.vecB)] = self._getDecagon(self.vecB, self.vecA)
        self.map_decagons[coords_tuple(self.vecC)] = self._getDecagon(self.vecC, self.vecA)
        self.map_decagons[coords_tuple(self.vecD)] = self._getDecagon(self.vecD, self.vecA)
        self.map_decagons[coords_tuple(self.vecE)] = self._getDecagon(self.vecE, self.vecA)
        self.map_decagons[coords_tuple(self.vecF)] = self._getDecagon(self.vecF, self.vecA)
        self.map_decagons[coords_tuple(-self.vecA)] = self._getDecagon(-self.vecA, -self.vecB)
        self.map_decagons[coords_tuple(-self.vecB)] = self._getDecagon(-self.vecB, -self.vecA)
        self.map_decagons[coords_tuple(-self.vecC)] = self._getDecagon(-self.vecC, -self.vecA)
        self.map_decagons[coords_tuple(-self.vecD)] = self._getDecagon(-self.vecD, -self.vecA)
        self.map_decagons[coords_tuple(-self.vecE)] = self._getDecagon(-self.vecE, -self.vecA)
        self.map_decagons[coords_tuple(-self.vecF)] = self._getDecagon(-self.vecF, -self.vecA)
        self.triangles = [self._getTriangle(self.vecA, self.vecB, self.vecC), \
                          self._getTriangle(self.vecA, self.vecC, self.vecD), \
                          self._getTriangle(self.vecA, self.vecD, self.vecE), \
                          self._getTriangle(self.vecA, self.vecE, self.vecF), \
                          self._getTriangle(self.vecA, self.vecF, self.vecB), \
                          self._getTriangle(self.vecC, self.vecB, -self.vecE), \
                          self._getTriangle(self.vecD, self.vecC, -self.vecF), \
                          self._getTriangle(self.vecE, self.vecD, -self.vecB), \
                          self._getTriangle(self.vecF, self.vecE, -self.vecC), \
                          self._getTriangle(self.vecB, self.vecF, -self.vecD), \
                          self._getTriangle(-self.vecA, -self.vecC, -self.vecB), \
                          self._getTriangle(-self.vecA, -self.vecD, -self.vecC), \
                          self._getTriangle(-self.vecA, -self.vecE, -self.vecD), \
                          self._getTriangle(-self.vecA, -self.vecF, -self.vecE), \
                          self._getTriangle(-self.vecA, -self.vecB, -self.vecF), \
                          self._getTriangle(-self.vecB, -self.vecC, self.vecE), \
                          self._getTriangle(-self.vecC, -self.vecD, self.vecF), \
                          self._getTriangle(-self.vecD, -self.vecE, self.vecB), \
                          self._getTriangle(-self.vecE, -self.vecF, self.vecC), \
                          self._getTriangle(-self.vecF, -self.vecB, self.vecD)]

    def _getDecagon(self, vec0, vec1):
        m10 = (vec0+vec1)/(1+self.__cosA) - vec0
        n10 = np.cross(vec0, m10)/math.sqrt(vec0.dot(vec0))
        decagon = []
        for i in range(0,10):
            v = vec0 + math.cos(math.pi/10+i*math.pi/5)*m10/math.cos(math.pi/10) + math.sin(math.pi/10+i*math.pi/5)*n10/math.cos(math.pi/10)
            key = coords_tuple(v)
            if key not in self.__joint:
                self.__joint[key] = v
            decagon.append(self.__joint[key])
        return decagon

    def _getTriangle(self, vec0, vec1, vec2):
        dec0 = self._getDecagon(vec0, vec1)
        dec1 = self._getDecagon(vec1, vec2)
        dec2 = self._getDecagon(vec2, vec0)
        v0 = self.__joint.get(coords_tuple(dec0[0]))
        v1 = self.__joint.get(coords_tuple(dec1[0]))
        v2 = self.__joint.get(coords_tuple(dec2[0]))
        if v0 is None or v1 is None or v2 is None:
            print 'error'
            return None
        return [v0, v1, v2]

def equals(a, b):
    if abs(a-b) < eps:
        return True
    else:
        return False

def writefacet(writefile, facet, scale):
    n = np.cross(facet[1]-facet[0],facet[2]-facet[0])
    n = n/math.sqrt(n.dot(n))
    for i in range(2,len(facet)):
        writefile.write('  facet normal %.10f %.10f %.10f\n' % (n[0], n[1], n[2]))
        writefile.write('    outer loop\n')
        writefile.write('      vertex %.10f %.10f %.10f\n' % (scale*facet[0][0], scale*facet[0][1], scale*facet[0][2]))
        writefile.write('      vertex %.10f %.10f %.10f\n' % (scale*facet[i-1][0], scale*facet[i-1][1], scale*facet[i-1][2]))
        writefile.write('      vertex %.10f %.10f %.10f\n' % (scale*facet[i][0], scale*facet[i][1], scale*facet[i][2]))
        writefile.write('    endloop\n')
        writefile.write('  endfacet\n')

#def less_than(vec1, vec2):
#    for i in range(0, len(vec1)):
#        if vec1[i] == vec2[i]: continue
#        if vec1[i] < vec2[i]: return True
#        if vec1[i] > vec2[i]: return False
#    return False
        
def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def coords_tuple(pt):
    return (round(pt[0],precision), round(pt[1],precision), round(pt[2],precision))

def udir(vec): # unit direction as a tuple
    tvec = np.asarray(vec)
    vlen = math.sqrt(tvec.dot(tvec))
    return (round(tvec[0]/vlen,precision), round(tvec[1]/vlen,precision), round(tvec[2]/vlen,precision))

def Read3DPenrose():
    r = t/2
    minx = float('inf')
    maxx = float('-inf')
    miny = float('inf')
    maxy = float('-inf')
    minz = float('inf')
    maxz = float('-inf')
    edges = []
    with open(input_filename, 'r') as readfile:
        for line in readfile:
            line = line.strip()
            ptInfo = line.split(',')
        
            x0 = float(ptInfo[0])
            y0 = float(ptInfo[1])
            z0 = float(ptInfo[2])
            x1 = float(ptInfo[3])
            y1 = float(ptInfo[4])
            z1 = float(ptInfo[5])
            
            minx = min([minx, x0, x1])
            miny = min([miny, y0, y1])
            minz = min([minz, z0, z1])
            maxx = max([maxx, x0, x1])
            maxy = max([maxy, y0, y1])
            maxz = max([maxz, z0, z1])
            
            edges.append([[x0,y0,z0],[x1,y1,z1]])

    trimmed_edges = [x for x in edges if not (x[0][0] < minx + trim or x[0][0] > maxx - trim or \
                                              x[0][1] < miny + trim or x[0][1] > maxy - trim or \
                                              x[0][2] < minz + trim or x[0][2] > maxz - trim or \
                                              x[1][0] < minx + trim or x[1][0] > maxx - trim or \
                                              x[1][1] < miny + trim or x[1][1] > maxy - trim or \
                                              x[1][2] < minz + trim or x[1][2] > maxz - trim)]
    map_vert = dict()
    map_vert_edges = dict()
    for i in range(0,len(trimmed_edges)):
        key0 = coords_tuple(trimmed_edges[i][0])
        key1 = coords_tuple(trimmed_edges[i][1])
        p0 = np.array(trimmed_edges[i][0])
        p1 = np.array(trimmed_edges[i][1])
        trimmed_edges[i] = [key0, key1]
        if key0 not in map_vert:
            map_vert[key0] = p0
            map_vert_edges[key0] = []
        if key1 not in map_vert:
            map_vert[key1] = p1
            map_vert_edges[key1] = []

        map_vert_edges[key0].append([key0, key1])
        map_vert_edges[key1].append([key1, key0])

    vecA = map_vert[trimmed_edges[0][1]]-map_vert[trimmed_edges[0][0]]
    vecB = None
    for edge in trimmed_edges:
        vecB = map_vert[edge[1]]-map_vert[edge[0]]
        if equals(vecA.dot(vecB), vecA.dot(vecA)/math.sqrt(5)):
            break

    if not (equals(vecA.dot(vecA), vecB.dot(vecB)) and equals(vecA.dot(vecB), vecA.dot(vecA)/math.sqrt(5))):
        print vecA.dot(vecA), vecB.dot(vecB), vecA.dot(vecB), vecA.dot(vecA)/math.sqrt(5)
        print 'error'
        return

    ico = Icosahedron(vecA, vecB, r)
#    print ico.verts()
    faces = []
    for edge in trimmed_edges:
        p0 = map_vert.get(edge[0])
        p1 = map_vert.get(edge[1])
        vec = ico.getDirection(p1-p0)
        if vec is None:
            print 'error'
            return
        dec = ico.getDecagon(vec)
        p_prev = dec[len(dec)-1]
        for i in range (0, len(dec)):
            faces.append([p0+dec[i],p1-2*vec+dec[i],p1-2*vec+p_prev,p0+p_prev])
            p_prev = dec[i]

    for key,pt in map_vert.iteritems():
        for tri in ico.triangles:
            faces.append([pt+tri[0], pt+tri[1], pt+tri[2]])
        adj_edges = map_vert_edges.get(key)
        empty_dirs = ico.verts()
        for edge in adj_edges:
            p0 = map_vert.get(edge[0])
            p1 = map_vert.get(edge[1])
            d = ico.getDirection(p1-p0)
            for i in range(0, len(empty_dirs)):  # instead of empty_dirs.remove(d)
                if (np.array_equal(empty_dirs[i], d)):
                    empty_dirs.pop(i)
                    break

        for d in empty_dirs:
            dec = list(ico.getDecagon(d))
            for i in range(0,10):
                dec[i] = pt + dec[i]
            faces.append(dec)


    with open(output_filename, 'w') as writefile:
        writefile.write('solid PENROSE3D\n')
        minv = [float('inf'), float('inf'), float('inf')]
        maxv = [float('-inf'), float('-inf'), float('-inf')]
        for f in faces:
            if len(f) < 3:
                print 'error'
                return
            writefacet(writefile, f, scale)
            for p in f:
                minv = np.minimum(minv, p)
                maxv = np.maximum(maxv, p)

        if plate_thickness > 0:
            offset = plate_thickness/2.0
            xmin = minv[0]
            ymin = minv[1]
            zmin = minv[2]
            xmax = maxv[0]
            ymax = maxv[1]
            zmax = maxv[2]
    ### bottom box
            writefacet(writefile, [np.array([xmin, ymin, zmin-offset]), np.array([xmin, ymax, zmin-offset]), np.array([xmax, ymax, zmin-offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmin-offset]), np.array([xmax, ymax, zmin-offset]), np.array([xmax, ymin, zmin-offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmin+offset]), np.array([xmax, ymax, zmin+offset]), np.array([xmin, ymax, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmin+offset]), np.array([xmax, ymin, zmin+offset]), np.array([xmax, ymax, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmin-offset]), np.array([xmax, ymin, zmin-offset]), np.array([xmax, ymin, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmin-offset]), np.array([xmax, ymin, zmin+offset]), np.array([xmin, ymin, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymin, zmin-offset]), np.array([xmax, ymax, zmin-offset]), np.array([xmax, ymax, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymin, zmin-offset]), np.array([xmax, ymax, zmin+offset]), np.array([xmax, ymin, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymax, zmin-offset]), np.array([xmin, ymax, zmin-offset]), np.array([xmin, ymax, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymax, zmin-offset]), np.array([xmin, ymax, zmin+offset]), np.array([xmax, ymax, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymax, zmin-offset]), np.array([xmin, ymin, zmin-offset]), np.array([xmin, ymin, zmin+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymax, zmin-offset]), np.array([xmin, ymin, zmin+offset]), np.array([xmin, ymax, zmin+offset])], scale)
    ### top box            
            writefacet(writefile, [np.array([xmin, ymin, zmax-offset]), np.array([xmin, ymax, zmax-offset]), np.array([xmax, ymax, zmax-offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmax-offset]), np.array([xmax, ymax, zmax-offset]), np.array([xmax, ymin, zmax-offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmax+offset]), np.array([xmax, ymax, zmax+offset]), np.array([xmin, ymax, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmax+offset]), np.array([xmax, ymin, zmax+offset]), np.array([xmax, ymax, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmax-offset]), np.array([xmax, ymin, zmax-offset]), np.array([xmax, ymin, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymin, zmax-offset]), np.array([xmax, ymin, zmax+offset]), np.array([xmin, ymin, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymin, zmax-offset]), np.array([xmax, ymax, zmax-offset]), np.array([xmax, ymax, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymin, zmax-offset]), np.array([xmax, ymax, zmax+offset]), np.array([xmax, ymin, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymax, zmax-offset]), np.array([xmin, ymax, zmax-offset]), np.array([xmin, ymax, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmax, ymax, zmax-offset]), np.array([xmin, ymax, zmax+offset]), np.array([xmax, ymax, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymax, zmax-offset]), np.array([xmin, ymin, zmax-offset]), np.array([xmin, ymin, zmax+offset])], scale)
            writefacet(writefile, [np.array([xmin, ymax, zmax-offset]), np.array([xmin, ymin, zmax+offset]), np.array([xmin, ymax, zmax+offset])], scale)
        
        writefile.write('endsolid PENROSE3D\n')



start_time = time.time()
Read3DPenrose()
print("--- %s seconds ---" % (time.time() - start_time))
