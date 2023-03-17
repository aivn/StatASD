#ifndef MAGNONS_HPP
#define MAGNONS_HPP

#include <aiwlib/sphere>
#include <aiwlib/gauss>
#include <aiwlib/mesh>

using namespace aiw;
//------------------------------------------------------------------------------
class MagnonsStochSrc{
	struct magnon_t {
		float phi0, A;  // начальная фаза и амплитуда
		Vecf<3> n, n_perp, k;
		Vecf<3> operator()(const Vecf<3> &r) const {
			return 	rotate(n_perp, n, sin(phi0+k*r))*A;
		}
	};
	Mesh<magnon_t, 3> magnons;
public:
	MagnonsStochSrc(int data_rank, int mg_split, float hTalpha){
		// magnons.init(Ind<3>(mg_split), Vecf<3>(), Vecf<3>(float(1<<data_rank)));
		magnons.init(Ind<3>(mg_split), Vecf<3>(), vecf(1<<data_rank));
		magnons.periodic = 7;

		sph_init_table(5);
		for(magnon_t &mg: magnons){
			mg.k[random()%3] = 2*M_PI;
			mg.n = sph_cell(rand()%sph_cells_num(5), 5);
			mg.n_perp = perp(mg.n);
			mg.phi0 = 2*M_PI*rand_alpha*rand();
			mg.A = sqrt(2*hTalpha)*rand_gauss();
			// mg.A = hTalpha;
		}
	}

	Vecf<3> operator ()(Vecf<3> r) const {
		r = magnons.coord2cell_bc(r);  Ind<3> pos = r;  r -= pos;  Vecf<3> H;
		// Vecf<4> w0[3]; for(int i=0; i<3; i++) w0[i] = interpolate_Bspline_weights<float>(r[i]);   
		// for(Ind<3> d: irange<3>(Ind<3>(4)))	H +=  w0[0][d[0]]*w0[1][d[1]]*w0[2][d[2]]*magnons[pos+d-ind(1)](r-d+ind(1));
		Vecf<2> w0[3]; for(int i=0; i<3; i++) w0[i] = interpolate_line_weights<float>(r[i]);   
		for(Ind<3> d: irange<3>(Ind<3>(2)))	H +=  w0[0][d[0]]*w0[1][d[1]]*w0[2][d[2]]*magnons[pos+d](r-d);
		return H;
	}
};
//------------------------------------------------------------------------------
class MagnonsStochSrcV2{
	Mesh<Vecf<8>, 3> magnons;
public:
	MagnonsStochSrcV2(int data_rank, int mg_split, float hTalpha){
		// magnons.init(Ind<3>(mg_split), Vecf<3>(), Vecf<3>(float(1<<data_rank)));
		magnons.init(Ind<3>(mg_split), Vecf<3>(), vecf(1<<data_rank));
		magnons.periodic = 7;
		magnons.interp = 0x111;
		
		sph_init_table(5);
		for(auto &mg: magnons){
			mg = Vecf<3>()|sph_cell(rand()%sph_cells_num(5), 5)|0|1; // M_PI*rand_alpha*rand()|1; // sqrt(2*hTalpha)*rand_gauss();
			mg[random()%3] = 2*M_PI*mg_split/(1<<data_rank);
		}
		magnons.dump(File("magnons.msh", "w"));
	}

	Vecf<3> operator ()(Vecf<3> r) const {
		Vecf<8> mg = magnons(r); Vecf<3> n = mg(3,4,5), n_perp = perp(n); n /= n.abs();
		return  rotate(n_perp, n, sin(mg(0,1,2)*r+mg[6]))*mg[7];
		// return  rotate(vecf(1,0,0), n, sin(mg(0,1,2)*r+mg[6]))*mg[7];
	}
};
//------------------------------------------------------------------------------
class MagnonsStochSrcV0{
	Vecf<3> k, n, n_perp; float phi0, A;   
public:
	MagnonsStochSrcV0(int data_rank, int, float hTalpha){
		for(int i=0; i<3; i++) if(rand()%2) k[i] = 2*M_PI/(1<<data_rank);
		A = sqrt(2*hTalpha)*rand_gauss();
		sph_init_table(5);
		n = sph_cell(rand()%sph_cells_num(5), 5);
		n_perp = perp(n);
		phi0 = rand()*rand_alpha*2*M_PI;
	}
	Vecf<3> operator ()(Vecf<3> r) const { return  rotate(n_perp, n, sin(phi0+k*r))*A; }
};
//------------------------------------------------------------------------------
#endif //MAGNONS_HPP
