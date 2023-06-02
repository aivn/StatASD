#ifndef MAGNONS_HPP
#define MAGNONS_HPP

#include <vector>
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
	MagnonsStochSrc(int data_rank, int mg_split, float hTalpha);
	
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
	MagnonsStochSrcV2(int data_rank, int mg_split, float hTalpha);
	
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
	MagnonsStochSrcV0(int data_rank, int mg_split, float hTalpha);
	
	Vecf<3> operator ()(Vecf<3> r) const { return  rotate(n_perp, n, phi0+k*r)*A; }
};
//------------------------------------------------------------------------------
class MagnonsStochSrcV3{
	Mesh<Vecf<3>, 3> magnons;
	Vecf<3> k, n, n_perp; float phi0, A;
public:
	MagnonsStochSrcV3(int data_rank, int mg_split, float hTalpha);
	
	Vecf<3> operator ()(Vecf<3> r) const {
		Vecf<3> mg = magnons(r);
		return  rotate(rotate(n_perp, n, phi0+k*r)*A, mg);
	}
};
//------------------------------------------------------------------------------
class MagnonsStochField{
	Mesh<Vecf<3>, 3> magnons;
	Vecf<3> k, n, n_perp; float phi0, A;
public:
	MagnonsStochField(int data_rank, Ind<3> mg_sz, float A_);
	
	Vecf<3> operator ()(Vecf<3> r) const {
		Vecf<3> mg = magnons(r);
		return  rotate(rotate(n_perp, n, phi0+k*r)*A, mg);
	}
};
//------------------------------------------------------------------------------
class MagnonsStochFieldGen{
	int data_rank;
	float halpha;
	
	struct node_t{
		Ind<3> k;
		float A, P; // амплитуда и интегральная вероятность
	};
	std::vector<node_t> table;
	float Pmax;
public:
	void init(int data_rank_, float T, float cL, float halpha_);
	float get_A(Ind<3> &k) const;
	MagnonsStochField operator ()() const;
};
//------------------------------------------------------------------------------
#endif //MAGNONS_HPP
