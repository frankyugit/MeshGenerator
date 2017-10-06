import json

class Mesh(object):

    def __init__(self, Nx, Ny, Nz, Lx, Ly, Lz):
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dx = Lx / Nx
        self.dy = Ly / Ny
        self.dz = Lz / Nz

        self.Ncell = Nx*Ny*Nz
        self.Nface = Nx*Ny*(Nz+1) + Ny*Nz*(Nx+1) + Nz*Nx*(Ny+1)
        self.Npoint = (Nx + 1) * (Ny + 1) * (Nz + 1)
        self.Ninnerface = 0

        self.Points = []
        self.genPoints()

        self.Faces = []
        self.genFaces()

        self.Cells = []
        self.genCells()

        self.FaceCells = []
        self.genFaceCells()




    def genPoints(self):
        for k in range(self.Nz + 1):
            for j in range(self.Ny + 1):
                for i in range(self.Nx + 1):
                    pCoordx = i * self.dx
                    pCoordy = j * self.dy
                    pCoordz = k * self.dz
                    Point = [pCoordx, pCoordy, pCoordz]
                    self.Points.append(Point)

    def genFaces(self):

        Nx = self.Nx
        Ny = self.Ny
        Nz = self.Nz

        # face generation along z direction
        for k in range(Nz + 1):
            for j in range(Ny):
                for i in range(Nx):
                    point1 = i+(Nx+1)*j+(Nx+1)*(Ny+1)*k
                    point2 = i+1+(Nx+1)*j+(Nx+1)*(Ny+1)*k
                    point3 = i+1+(Nx+1)*(j+1)+(Nx+1)*(Ny+1)*k
                    point4 = i+(Nx+1)*(j+1)+(Nx+1)*(Ny+1)*k
                    face = [point1,point2,point3,point4]
                    self.Faces.append(face)


        # face generation along x direction
        for i in range(Nx + 1):
            for k in range(Nz):
                for j in range(Ny):
                    point1 = i+(Nx+1)*j+(Nx+1)*(Ny+1)*k
                    point2 = i+(Nx+1)*(j+1)+(Nx+1)*(Ny+1)*k
                    point3 = i+(Nx+1)*(j+1)+(Nx+1)*(Ny+1)*(k+1)
                    point4 = i+(Nx+1)*j+(Nx+1)*(Ny+1)*(k+1)
                    face = [point1,point2,point3,point4]
                    self.Faces.append(face)

        # face generation along y direction
        for j in range(Ny + 1):
            for i in range(Nx):
                for k in range(Nz):
                    point1 = i+(Nx+1)*j+(Nx+1)*(Ny+1)*k
                    point2 = i+(Nx+1)*j+(Nx+1)*(Ny+1)*(k+1)
                    point3 = i+1+(Nx+1)*j+(Nx+1)*(Ny+1)*(k+1)
                    point4 = i+1+(Nx+1)*j+(Nx+1)*(Ny+1)*k
                    face = [point1,point2,point3,point4]
                    self.Faces.append(face)

    def genCells(self):
        Nx = self.Nx
        Ny = self.Ny
        Nz = self.Nz
        for k in range(Nz):
            for j in range(Ny):
                for i in range(Nx):
                    face1 = i + Nx*j + Nx*Ny*k
                    face2 = i + Nx*j + Nx*Ny*(k+1)
                    face3 = j + Ny*k + Ny*Nz*i + Nx*Ny*(Nz+1)
                    face4 = j + Ny*k + Ny*Nz*(i+1) + Nx*Ny*(Nz+1)
                    face5 = k + Nz*i + Nz*Nx*j + Nx*Ny*(Nz+1) + Ny*Nz*(Nx+1)
                    face6 = k + Nz*i + Nz*Nx*(j+1) + Nx*Ny*(Nz+1) + Ny*Nz*(Nx+1)
                    cell = [face1, face2, face3, face4, face5, face6]
                    self.Cells.append(cell)

    def genFaceCells(self):
        Nface = self.Nface
        FaceCells = []
        for i in range(Nface):
            FaceCells.append([-1,-1])
        celli = 0
        for cell in self.Cells:
            for face in cell:
                facecell = FaceCells[face]
                if facecell[0] == -1:
                    facecell[0] = celli
                elif facecell[1] == -1:
                    facecell[1] = celli
                else:
                    print("Face %d got more 2 cell neighbours" % face)
            celli = celli + 1
        self.FaceCells = FaceCells


    def writeMesh(self, filename, pathname='./'):
        print("Mesh Writing...")
        fid = open(pathname + filename, 'w')
        fid.write("Mesh\n")
        fid.write("{\n")
        fid.write("Ncell:%d\nNface:%d\nNpoint:%d\nNinnerface:%d\n" % (self.Ncell, self.Nface, self.Npoint, self.Ninnerface))
        fid.write("}\n")

        fid.write("Points\n")
        fid.write("{\n")
        for point in self.Points:
            fid.write("(%f, %f, %f)\n" % (point[0], point[1], point[2]))
        fid.write("}\n")

        fid.write("Faces\n")
        fid.write("{\n")
        for face in self.Faces:
            fid.write("(%d, %d, %d, %d)\n" % (face[0], face[1], face[2], face[3]))
        fid.write("}\n")

        fid.write("Cells\n")
        fid.write("{\n")
        for cell in self.Cells:
            fid.write("(%d, %d, %d, %d, %d, %d)\n" % (cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]))
        fid.write("}\n")

        fid.write("FaceCells\n")
        fid.write("\n")
        fid.write("{\n")
        for facecell in self.FaceCells:
            fid.write("(%d, %d)\n" % (facecell[0], facecell[1]))
        fid.write("}\n")

        fid.close()
        print("Done")

    def writeJSONMesh(self, filename, pathname='./'):
        print("Json format Mesh Writing...")

        Jmesh = {
            'Mesh': {'Ncell': self.Ncell, 'Nface': self.Nface, 'Npoint': self.Npoint, 'Ninnerface': self.Ninnerface},
            'Points': self.Points,
            'Faces': self.Faces,
            'Cells': self.Cells,
            'FaceCells': self.FaceCells

        }

        with open(filename, 'w') as f:
            json.dump(Jmesh, f)
        print("Done")






if __name__=="__main__":
    import sys

    Nx = 5
    Ny = 5
    Nz = 1

    Lx = 10
    Ly = 10
    Lz = 0.1

    mesh = Mesh(Nx, Ny, Nz, Lx, Ly, Lz)

    mesh.writeMesh("mesh")
    mesh.writeJSONMesh("mesh.json")
    sys.exit()