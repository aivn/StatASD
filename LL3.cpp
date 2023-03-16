#define LL3_CPP
#include <omp.h>
#include "LL3.hpp"
#include "magnons.hpp"

#ifdef CCC
const char *mode = "CCC";
#endif
#ifdef VCC
const char *mode = "VCC";
#endif
#ifdef FCC
const char *mode = "FCC";
#endif

const double _3 = 1./3, _6 = 1./6;
//------------------------------------------------------------------------------
void Model::init(const char *path_){
	path = path_; M0 /= M0.abs();

	sph_init_table(ind(5, f_rank).max());
	if(f_use){ f.init(f_rank, 0, 1); f_eq.init(f_rank, 0, 1); f_eq.fill(0.f); }

	for(int k=0; k<4; k++) data[k].init(data_rank);
	if(stoch0){ for(size_t i=0; i<data[0].size(); i++) for(int k=0; k<cell_sz; k++){
				data[0][i].m[k] = sph_cell(rand()%sph_cells_num(5), 5)+M0;  data[0][i].m[k] /= data[0][i].m[k].abs();
			}
	} else if(helic0) {
		for(size_t i=0; i<data[0].size(); i++){
			double phi = 2*M_PI*data[0].offset2pos(i)[0]/data[0].bbox()[0];
			Vecf<3> m0; m0[0] = cos(phi); m0[1] = sin(phi); m0[2] = helic_z;
			data[0][i].init(m0/m0.abs());
		}
	} else for(size_t i=0; i<data[0].size(); i++) data[0][i].init(M0);
	t = 0.;
	ftm = File("%tm0.dat", "w", path_); ftm("#:t mx my mz Hx Hy Hz\n% % %\n", t, data[0][ind(0,0,0)].m[0], Hexch(0, 0, data[0].get_nb(0, 7), 0));
	ftvals = File("%tvals.dat", "w", path_);
	ftvals("#:t M Mx My Mz M2x M2y M2z W Wexch Wext Wanis c20 c21 c22  Q1 Q2 Q3  eta eta2 eta3 eta4\n");

	// M = M2 = vec(0., 0., 1.);
	// W[1] = -nb_sz/2*J; W[2] = -Hext[2]; W[3] = -K*nK[2]*nK[2]; W[0] = W[1]+W[2]+W[3];
	// WOUT(data_rank, cell_sz, nb_sz, mode);
	calc_av();

	if(calc_cl) chain_lambda.resize(1<<(data_rank-1));
}
//------------------------------------------------------------------------------
void Model::calc_av(){  // считаем средние значения W, M, M2
	size_t sz = data[0].size();
	double Mx = 0, My = 0, Mz = 0, M2x = 0, M2y = 0, M2z = 0, Wx = 0, We = 0, Wa = 0, C20 = 0, C21 = 0, C22 = 0, Q20 = 0, Q21 = 0, Q22 = 0, eta1 = 0, eta2 = 0, eta3 = 0, eta4 = 0; 
#pragma omp parallel for reduction(+:Mx,My,Mz,M2x,M2y,M2z,Wx,We,Wa,C20,C21,C22,Q20,Q21,Q22,eta1,eta2,eta3,eta4) if(parallel)
	for(size_t i=0; i<sz; ++i){
		ZCubeNb<3> nb = data[0].get_nb(i, 7);
		for(int k=0; k<cell_sz; k++){
			const Vecf<3> &m0 = data[0][i].m[k];
			Mx += m0[0];        My += m0[1];        Mz += m0[2];
			M2x += m0[0]*m0[0]; M2y += m0[1]*m0[1]; M2z += m0[2]*m0[2];
			float e = Hexch(0, i, nb, k)*m0;  Wx -= e*.5f;  e /= nb_sz;
			eta1 += e;  eta2 += e*e;  eta3 += e*e*e;  eta4 += e*e*e*e;
			We -= Hext*m0;
			float nKm = nK*m0; Wa -= K*nKm*nKm;

			Vecf<corr2_types> C2, Q2;
			for(int ci=0; ci<corr2_types; ci++)
				for(int j=0; j<corr2_sz[ci]; j++){
					const Ind<2> &link = corr2_table[ci*cell_sz+k][j];
					const Link &l0 = nb_pos[k][link[0]], &l1 = nb_pos[k][link[1]];
					const auto &mi = data[0][i+nb[l0.off]].m[l0.cell], &mk = data[0][i+nb[l1.off]].m[l1.cell];
					C2[ci] += mi*mk;
					Q2[ci] += mi*(m0%(m0%mk));
				}
			C20 += C2[0]; C21 += C2[1]; Q20 += Q2[0]; Q21 += Q2[1]; if(corr2_types>=3){ C22 += C2[2]; Q22 += Q2[2]; }
		}
	}

	M  = vec(Mx, My, Mz)/(sz*cell_sz);
	M2 = vec(M2x, M2y, M2z)/(sz*cell_sz);
	W  = vec(Wx+We+Wa, Wx, We, Wa)/(sz*cell_sz);
	corr2[0] = C20/(sz*cell_sz*corr2_sz[0]);
	corr2[1] = C21/(sz*cell_sz*corr2_sz[1]);
	if(corr2_types>=3) corr2[2] = C22/(sz*cell_sz*corr2_sz[2]);
	Q[0] = Q20/(sz*cell_sz*corr2_sz[0]);
	Q[1] = Q21/(sz*cell_sz*corr2_sz[1]);
	if(corr2_types>=3) Q[2] = Q22/(sz*cell_sz*corr2_sz[2]);
	eta = vecf(eta1, eta2, eta3, eta4)/(sz*cell_sz);
}
//------------------------------------------------------------------------------
void Model::calc(int steps){
	float sghT = sqrt(2*dt*alpha*T); // 2*sqrt(dt*alpha*T);
	size_t sz = data[0].size();
	for(int Nt=0; Nt<steps; Nt++){
#pragma omp parallel for if(parallel)
		for(size_t i=0; i<sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(0, i, nb, k);
				WOUT(i, k, Hex); 
				Vecf<3> H = Hex + Hext + 2*K*nK*m0*nK; 
				Vecf<3> dmdt = m0%(-gamma*H -alpha*m0%H);
				m1 = m0 + .5f*dt*dmdt;
				dm = .5f*dmdt;
			}
			// exit(0);
		}
		//-----------
#pragma omp parallel for if(parallel)
		for(size_t i=0; i<sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &m2 = data[2][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(1, i, nb, k);
				Vecf<3> H = Hex + Hext + 2*K*nK*m1*nK; 
				Vecf<3> dmdt = m1%(-gamma*H -alpha*m1%H);
				m2 = m0 + .5f*dt*dmdt;
				dm += dmdt;
			}
		}
		//------------
#pragma omp parallel for if(parallel)
		for(size_t i=0; i<sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &m2 = data[2][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(2, i, nb, k);
				Vecf<3> H = Hex + Hext + 2*K*nK*m2*nK; 
				Vecf<3> dmdt = m2%(-gamma*H -alpha*m2%H);
				m1 = m0 + dt*dmdt;
				dm += dmdt;
			}
		}
		//-------------
		// делаем магнон
		/*
		Vecf<3> km(random()%nm_max, random()%nm_max, random()%nm_max); km *= 2*M_PI/(1<<data_rank);
		float phi0 = rand_alpha*random()*2*M_PI;
		Vecf<3> nm = sph_cell(rand()%sph_cells_num(5), 5), nm_perp = perp(nm);
		*/
		MagnonsStochSrcV2 mg(data_rank, mg_split, dt*T*alpha);

		// unsigned int seed = 0;		
		// #pragma omp parallel for firstprivate(seed) if(parallel)
#pragma omp parallel for if(parallel)
		for(size_t i=0; i<sz; ++i){
			// rand_init(seed, omp_get_thread_num());
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(1, i, nb, k);
				Vecf<3> H = Hex + Hext + 2*K*nK*m1*nK; 
				Vecf<3> dmdt = m1%(-gamma*H -alpha*m1%H);
				m0 += dt*_3*(dm + .5*dmdt);
				// m0 = gauss_rotate(m0/m0.abs(), sghT, seed); // old variant, errro
				m0 /= m0.abs(); // data[2][i].m[k] = m0;

				// Vecf<3> Hm = rotate(nm_perp, nm, sin(phi0+km*coord[k]+km*zoff2pos<3>(i, data_rank)));
				
				// m0 = rotate(m0, rand_gaussV<3, float>(seed)*sghT);
				// m0 = rotate(m0, Hm*sghT);
				m0 = rotate(m0, mg(zoff2pos<3>(i, data_rank)+coord[k]));
			}
		}
		//---------------
		/*
		unsigned int seed = 0;	
#pragma omp parallel for firstprivate(seed) if(parallel)
		for(size_t i=0; i<sz; ++i){
			rand_init(seed, omp_get_thread_num());
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m2 = data[2][i].m[k];
				Vecf<3> Hex = Hexch(2, i, nb, k);
				Vecf<3> H = Hex + Hext + 2*K*nK*m2*nK;
				// m0 = rotate(m0, rand_gaussV<3, float>(seed)*sghT*sqrtf(1.f-expf((H*m0+mu)/T)));
				m0 = rotate(m0, rand_gaussV<3, float>(seed)*sqrtf(2*dt*alpha*(T+A*(8-H*m0))));
			}
		}
		*/
		
		t += dt; ftm("% % %\n", t, data[0][0].m[0], Hexch(0, 0, data[0].get_nb(0, 7), 0));
		if(calc_eq || Nt==steps-1) calc_av();
		if(calc_eq){ Meq += M; Mabs_eq += M.abs(); M2eq += M2; Weq += W; corr2eq += corr2; Qeq += Q; eta_eq += eta; eq_count++; }
		if(calc_eq && calc_cl) calc_chain_lambda();
	}
	ftvals("% %   %       %       %      %      %     %\n", t, M.abs(), M, M2, W, corr2, Q, eta).flush();
}
//------------------------------------------------------------------------------
void Model::calc_chain_lambda(){
	int thN = omp_get_max_threads();
	// printf("thN=%i\n", thN);
	std::vector<Vec<4> > tmp_cl[thN]; for(int i=0; i<thN; i++) tmp_cl[i].resize(chain_lambda.size());
	size_t sz = data[0].size(); int  zsz =  1<<data_rank, hsz = zsz/2;
#pragma omp parallel for if(parallel)
	for(size_t i=0; i<sz; i++){ // цикл по XYZ
		// printf("thN=%i th=%i\n", thN, omp_get_thread_num());
		std::vector<Vec<4> > &cl = tmp_cl[omp_get_thread_num()];
		const Vecf<3> &mi = data[0][i].m[0]; Ind<3> pos = data[0].offset2pos(i); int z0 = pos[2];
		for(int j=0; j<hsz; j++){
			pos[2] = (z0+j+1)%zsz; float eta = mi*data[0][data[0].pos2offset(pos)].m[0], eta_k = eta;
			for(int k=0; k<4; k++){ cl[j][k] += eta_k; eta_k *= eta; } 
		} 
	} // конец цикла по XYZ
	double div_ = 1./(1<<3*data_rank);
	for(int i=0; i<thN; i++) for(int j=0, sz=chain_lambda.size(); j<sz; j++) chain_lambda[j] += div_*tmp_cl[i][j];
}
//------------------------------------------------------------------------------
void Model::clean_av_eq(){
	eq_count = 0; Meq = vec(0.); Mabs_eq = 0.; M2eq = vec(0.); Weq = vec(0.); corr2eq = vec(0.); Qeq = vec(0.); eta_eq = vec(0.);
}
void Model::calc_f(){
	if(!f_use) return;
	size_t sz = data[0].size();
	f.fill(0.f); 
	double  df = 1/(4*M_PI/f.size()*sz*cell_sz);
	
	for(size_t i=0; i<sz; ++i) for(int k=0; k<cell_sz; k++) f[Vec<3>(data[0][i].m[k])] += df;
	if(calc_eq){
		for(int i=0; i<int(f.size()); i++) f_eq[i] += f[i]; 
		eq_f_count++;
	}
}
//------------------------------------------------------------------------------
void Model::finish(){
	if(eq_count){
		Weq /= eq_count;
		Meq /= eq_count;
		M2eq /= eq_count;
		Mabs_eq /= eq_count;
		corr2eq /= eq_count;
		Qeq /= eq_count;
		eta_eq /= eq_count;
		if(f_use) for(int i=0; i<int(f_eq.size()); i++) f_eq[i] /= eq_f_count;
		if(calc_cl){
			File fcl("%cl.dat", "w", path.c_str()); fcl.printf("#:n eta1 eta2 eta3 eta4\n");
			for(int i=0, sz = chain_lambda.size(); i<sz; i++) fcl("% %\n", i, chain_lambda[i]/eq_count);
		}
	} else {
		Weq = W;
		Meq = M;
		M2eq = M2;
		Mabs_eq = M.abs();
		if(f_use) f_eq = f.copy();
		corr2eq = corr2;
		Qeq = Q;
		eta_eq = eta;
	}
	if(f_use) f_eq.dump(File("%f_eq.sph", "w", path.c_str()));
}
//------------------------------------------------------------------------------
// void Model::dump(aiw::IOstream &S);
// void Model::load(aiw::IOstream &S);
//------------------------------------------------------------------------------
