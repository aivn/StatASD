#include <omp.h>
#include "SCL.hpp"

using namespace aiw;
//------------------------------------------------------------------------------
void Model::init(const char *path){
	for(int i=0; i<4; i++) data[i].resize(mesh_sz.prod(), Vecf<3>(0.f, 0.f, 1.f));
	links.resize(data[0].size()*3);
		
	calc_off = CalcOff<3>(Vec<3, int64_t>(mesh_sz), periodic);
	Nth = omp_get_max_threads();
	randN01.init(Nth);

	std::uniform_real_distribution<double> u(0., 1.);
	std::mt19937 gen; std::random_device rd;  gen.seed(rd());
	for(size_t i=0, sz=links.size(); i<sz; i++) links[i] = u(gen)<sparse;
	
	ftvals.open(std::string(path)+"tvals.dat");
	W[1] = 0; for(int64_t off=0, sz=data[0].size(); off<sz; off++) W[1] -= Hexch(0, off)[2];
	W[1] /= 2*data[0].size(); W[2] = -Hext[2]; W[3] = -K*nK[2]*nK[2]; W[0] = W[1] + W[2] + W[3];
	n_b = -W[1]*2;
	ftvals<<"#:t M Mx My Mz M2x M2y M2z eta W Wexch Wext Wanis\n0 1 0 0 1 0 0 1 1 "<<W<<'\n';
}
//------------------------------------------------------------------------------
const float _3 = 1./3;
void Model::calc(int steps){
	int64_t data_sz = data[0].size();
	float Kx2 = 2*K,  sghT = sqrt(2*dt*alpha);
	for(int Nt=0; Nt<steps; Nt++){
#pragma omp parallel for 
		for(int64_t i=0; i<data_sz; ++i){
			Vecf<3> &m0 = data[0][i], &m1 = data[1][i], &dm = data[3][i];
			Vecf<3> Hex = Hexch(0, i);
			Vecf<3> H = Hex + Hext + nK*m0*Kx2*nK; 
			Vecf<3> dmdt = m0%(-gamma*H -alpha*m0%H);
			m1 = m0 + .5f*dt*dmdt;
			dm = .5f*dmdt;
		}
		//-----------
#pragma omp parallel for 
		for(int64_t i=0; i<data_sz; ++i){
			Vecf<3> &m0 = data[0][i], &m1 = data[1][i], &m2 = data[2][i], &dm = data[3][i];
			Vecf<3> Hex = Hexch(1, i);
			Vecf<3> H = Hex + Hext + nK*m1*Kx2*nK; 
			Vecf<3> dmdt = m1%(-gamma*H -alpha*m1%H);
			m2 = m0 + .5f*dt*dmdt;
			dm += dmdt;
		}
		//------------
#pragma omp parallel for
		for(int64_t i=0; i<data_sz; ++i){
			Vecf<3> &m0 = data[0][i], &m1 = data[1][i], &m2 = data[2][i], &dm = data[3][i];
			Vecf<3> Hex = Hexch(2, i);
			Vecf<3> H = Hex + Hext + nK*m2*Kx2*nK; 
			Vecf<3> dmdt = m2%(-gamma*H -alpha*m2%H);
			m1 = m0 + dt*dmdt;
			dm += dmdt;
		}
		//-------------
#pragma omp parallel for 
		for(int64_t i=0; i<data_sz; ++i){
			int thID = omp_get_thread_num();
			Vecf<3> &m0 = data[0][i], &m1 = data[1][i], &dm = data[3][i];
			Vecf<3> Hex = Hexch(1, i);
			Vecf<3> H = Hex + Hext + nK*m1*Kx2*nK; 
			Vecf<3> dmdt = m1%(-gamma*H -alpha*m1%H);
			m0 += dt*_3*(dm + .5f*dmdt);
			m0 /= m0.abs();

			m0 = rotate(m0, randN01.V<3>(thID)*sghT);  // ??? поворот или приращение ???
		}
		t += dt; 
		if(calc_eq || Nt==steps-1){
			double Mx = 0, My = 0, Mz = 0, M2x = 0, M2y = 0, M2z = 0, Wexch = 0, Wan = 0;
#pragma omp parallel for reduction(+:Mx, My, Mz, M2x, M2y, M2z, Wexch, Wan)			
			for(int64_t i=0; i<data_sz; ++i){
				Vecf<3> &m0 = data[0][i], Hex = Hexch(0, i);
				Mx += m0[0];  M2x += m0[0]*m0[0]; 
				My += m0[1];  M2y += m0[1]*m0[1]; 
				Mz += m0[2];  M2z += m0[2]*m0[2];
				float MnK = m0*nK;
				Wexch -= m0*Hex*.5f; Wan -= K*MnK*MnK;
			}
			M = vec(Mx, My, Mz)/data_sz;
			M2 = vec(M2x, M2y, M2z)/data_sz;
			W[1] = Wexch/data_sz; W[2] = -M*Hext; W[3] = Wan/data_sz; W[0] = W[1]+W[2]+W[3];
			eta = -W[1]/n_b*2;
			if(calc_eq){
				Meq += M; Mabs_eq += M.abs();
				M2eq += M2; Weq += W; eta_eq += eta;
				eq_count++;
			}
		}
	}
	ftvals<<t<<' '<<M.abs()<<' '<<M<<' '<<M2<<' '<<eta<<' '<<W<<std::endl;
}
//------------------------------------------------------------------------------
void Model::finish(){
	if(eq_count){
		Weq /= eq_count;
		Meq /= eq_count;
		M2eq /= eq_count;
		Mabs_eq /= eq_count;
		eta_eq /= eq_count;
	} else {
		Weq = W;
		Meq = M;
		M2eq = M2;
		Mabs_eq = M.abs();
		eta_eq = eta;
	}
}
//------------------------------------------------------------------------------
