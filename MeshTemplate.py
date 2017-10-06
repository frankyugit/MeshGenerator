from string import Template

mesh_template = Template(
"""
mesh
{
Ncell:$Ncell;Nface:$Nface;Npoint:%Npoint;Ninnerface:%Ninnerface
}
points
{
%points
}
faces
{
()
}
Cell
{

}
faceCell
{
()
}
"""
)