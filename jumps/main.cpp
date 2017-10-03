#include <dolfin.h>
#include <ufc.h>                                                                                                                                    
#include <dolfin/fem/UFC.h> 
//#include <dolfin/mesh/UnitSquare.h>
//#include <NSEResMomentum2D.h>
#include <basic2D.h>
using namespace dolfin;
int main()
{
//UnitSquare sq(2,2);
Mesh sq("us.xml");
Function u;
Function u2;
PETScVector uX;
PETScVector uX2;
basic2DLinearForm *form = new basic2DLinearForm(u);
File f("unitsq.pvd");
u.init(sq, uX, *form, 1);
u2.init(sq, uX2, *form, 1);

MeshFunction<real> meshfun; 
meshfun.init(sq, 2); 



MeshFunction<dolfin::uint> partitions;
partitions.init(sq,sq.topology().dim());
for(CellIterator c(sq); !c.end(); ++c)
	{
	int num = 0;
	if (c->midpoint().x() > 0.5)
		num = num +1;
	if (c->midpoint().y() > 0.5)
		num = num +2;
	
	partitions.set(c->index(), num);
}

if (dolfin::MPI::numProcesses() > 1)
{
sq.distribute(partitions);
sq.renumber();
}


UFC ufc(form->form(),sq, form->dofMaps());
Cell c(sq,0);
uint local_dim = c.numEntities(0);
uint *idx = new uint[local_dim];
uint *id = new uint[local_dim];
real *block = new real[local_dim];



for(CellIterator c(sq); !c.end(); ++c)
{
	uint ii = 0;
	uint jj = 0;


 ufc.update(*c, sq.distdata());
 (form->dofMaps())[1].tabulate_dofs(idx, ufc.cell, c->index()); 
// for(VertexIterator v(*c); !v.end(); ++v, ii++)
	 {
//		if (!sq.distdata().is_ghost(v->index(), 0))
		{
			block[jj] = c->midpoint().x() > 0.5;
			id[jj++]= idx[ii];      
			uX.set(block,jj,id);
		}
	}
                                                  
	
//	std::cout << MPI::processNumber() << " : " << sq.distdata().get_cell_global(c->index()) << std::endl;
	} 

uX.apply();
delete[] idx;
delete[] id;
delete[] block;

f << u;




//uX.disp();

Assembler a(sq);
a.assemble(uX2,*form, true);

uX2.disp();



 local_dim = c.numEntities(0);
 idx = new uint[local_dim];
 id = new uint[local_dim];
 block = new real[local_dim];
 


for(CellIterator c(sq); !c.end(); ++c)
{
        uint ii = 0;
        uint jj = 0;

        ufc.update(*c, sq.distdata());
        (form->dofMaps())[1].tabulate_dofs(idx, ufc.cell, c->index());
        uX2.get(block,1,idx);
        meshfun.set(*c,block[jj]);
}

delete[] idx;
delete[] id;
delete[] block;




File ff("uX2.pvd");
ff << meshfun;
//std::cout << uX.size()<< std::endl;


//if (dolfin::MPI::numProcesses() > 1)
MPI_Finalize();
}
