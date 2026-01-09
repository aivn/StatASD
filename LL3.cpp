// #define LL3_CPP
#include <fstream>
#include <map>
#include <omp.h>
#include <aiwlib/llbe>
#include "LL3.hpp"
// #include "Z2.hpp"

#ifdef MAGNONS
#include "magnons.hpp"
#endif // MAGNONS

#ifdef SC
const char *mode = "SC";
#endif
#ifdef BCC
const char *mode = "BCC";
#endif
#ifdef FCC
const char *mode = "FCC";
#endif

const double _3 = 1./3, _6 = 1./6;

inline double pow2(double x){ return x*x; }
float calc_Mu(float M){
	float M2 = M*M, tmp = 0.775383f + (0.109185f + (0.289114f + (-0.871214f + (1.85968f + (-1.71306f + 0.550912f*M2 )*M2 )*M2 )*M2 )*M2 )*M2;
	return (1-M2*tmp*tmp)/3.;
}
inline float Upsilon(float M, float eta){
	if(fabs(M)>.9999999f) return 0;
	float M2 = M*M, zeta = (eta-M2)/(1-M2), zeta2 = zeta*zeta, Mz3 = M+zeta; Mz3 *= Mz3*Mz3;
	float delta = 0.1223f*Mz3*Mz3*(1.f-1.6f*M-2.49f*zeta+0.71f*M2+1.71f*zeta2+2.47f*M*zeta-0.12f*M2*M-0.21f*zeta2*zeta-0.34f*M2*zeta-1.17f*M*zeta2);
	return 3*calc_Mu(M)*calc_Mu(zeta)/(1+zeta)*(1+delta);
}
/*
inline double Upsilon(double M, double eta){
	M = fabs(M); double M2 = M*M, zeta = M2<1? (eta-M2)/(1-M2): 0, zeta2 = zeta*zeta;
	double mu_p = (1-pow2( M*( 0.775383 + (0.109185 + (0.289114 + (-0.871214 + (1.85968 + (-1.71306 + 0.550912*M2 )*M2 )*M2 )*M2 )*M2 )*M2 ) ))/3.;
    return (1-zeta)*mu_p + 0.0467246*(1-M)*zeta*(1-zeta)*
		(1 - 7.96819*M - 0.939775*zeta + 14.0242*M*zeta + 20.9323*M2 
		 + 7.65054*zeta2 - 14.5937*M2*zeta - 5.35767*M*zeta2 
		 - 13.542*M2*M - 3.72091*zeta2*zeta);
}
*/
//------------------------------------------------------------------------------
void Model::init(const char *path_){
	double t0 = omp_get_wtime();
	if(threads!=0) omp_set_num_threads(threads);
	Nth = omp_get_max_threads();
	randN01.init(Nth);
	path = path_; // M0 /= M0.abs();  

	sph_init_table(ind(5, f_rank).max());
	if(f_rank>=0){ f.init(f_rank, 0, 1); f_eq.init(f_rank, 0, 1); f_eq.fill(0.f); f_buf.resize(Nth*f.size()); }
	if(fz_sz){ fz.resize(fz_sz); fz_eq.resize(fz_sz, 0.f); fz_buf.resize(Nth*fz_sz); }
	if(f_eta_sz){ f_eta.resize(f_eta_sz); f_eta_eq.resize(f_eta_sz, 0.f); f_eta_buf.resize(Nth*f_eta_sz); }

	t = 0.;   T_sc = T; Tsc_eq = 0;

	if(Ms_start>=0 && Ms_start<data_rank){
		Ms_arr.resize(data_rank-1);
		for(int i=0; i<data_rank-1; i++) Ms_arr[i].init(Ind<3>(1<<(data_rank-i)), Vec<3>(), Vec<3>(double(1<<data_rank)));
	} else { Ms_arr.resize(1); Ms_arr[0].init(Ind<3>(1)); }
	
	if(corr_max){
		if(corr_max<0) corr_max = 1<<(data_rank-1);
		corr.resize(corr_max); corr_eq.resize(corr_max); corr_buf.resize(Nth*corr_max);  corr_direct_offs.resize(size_t(corr_max)*corr_direct_sz*(1<<data_rank*3));
		for(size_t i=0, sz=corr_direct_offs.size(); i<sz; i++){  // usage:   corr_direct_offs[(cID*corr_max + l)*corr_direct_sz + d]
			int d = i%corr_direct_sz, r = 1+i/corr_direct_sz%corr_max, f = i/(corr_direct_sz*corr_max);
			Ind<3> pos = data[0].offset2pos(f)+corr_direct[d]*r;  for(int k=0; k<3; k++) pos[k] = (pos[k]+(1<<data_rank))%(1<<data_rank);
			corr_direct_offs[i] = data[0].pos2offset(pos);
		}
	}

	// rand_init();
#ifdef CALC_Q
	Q_buf.resize(Nth*Q_sz);  eta_k_buf.resize(Nth*Q_sz*4);
#endif  // CALC_Q
	
#ifdef MAGNONS
	MG.init(data_rank, T, cL, dt*alpha);
#endif // MAGNONS


	fftw.init(Ind<3>(1<<data_rank), 3);
	fspectrum.open(std::string(path)+"spectrum.dat"); fspectrum<<"#:t nu ampl\n";

	rt_init += omp_get_wtime()-t0;
}
//------------------------------------------------------------------------------
void Model::start_gauss(){
	double t0 = omp_get_wtime();
	for(int k=0; k<4; k++) data[k].init(data_rank);
	Natoms = data[0].size()*cell_sz;
	float M0abs = M0.abs();
	if(M0abs>=.99){
		M0 /= M0abs;
		for(int i=0, sz=data[0].size(); i<sz; i++) data[0][i].init(M0);
	} else {
		Vecf<3> p; if(M0abs>1e-3f) p = LLBE::invL(M0abs)/M0abs*M0;
		std::map<float, int> ftable; float s = 0; for(int i=0, sz=sph_vertex_num(5); i<sz; i++){ s += expf(sph_vert(i, 5)*p); ftable[s] = i; }
		std::random_device rd; std::mt19937 gen;  gen.seed(rd());  std::uniform_real_distribution<float>  rnd(0, s);
		for(size_t i=0; i<data[0].size(); i++) for(int k=0; k<cell_sz; k++) data[0][i].m[k] = sph_vert(ftable.lower_bound(rnd(gen))->second, 5);
	}
	rt_init += omp_get_wtime()-t0;
	init_conditions = true;  // флаг задания н.у.
	open_tvals();
}
void Model::start_helic(){
	double t0 = omp_get_wtime();
	float M0abs = M0.abs(); if(M0abs>=1 || M0abs<1e-3){  start_gauss();  return; }
	Vecf<3> n = M0/M0abs, m0 = M0 + perp(n)*sqrtf(1-M0*M0);
	for(int k=0; k<4; k++) data[k].init(data_rank);
	Natoms = data[0].size()*cell_sz;
	for(size_t i=0, sz=data[0].size(); i<sz; i++){
		float z = data[0].offset2pos(i)[2];
		for(int k=0; k<cell_sz; k++) data[0][i].m[k] = rotate(m0, n, 2*M_PI*(z+coord[k][2])*helic_n/data[0].bbox()[2]);
	}
	rt_init += omp_get_wtime()-t0;
	init_conditions = true;  // флаг задания н.у.
	open_tvals();
}
//------------------------------------------------------------------------------
void Model::open_tvals(){
	t = 0;
	if(out_tm0){
		ftm.close(); ftm = File("%tm0.dat", "w", path);
		ftm("#:t mx my mz Hx Hy Hz\n% % %\n", t, data[0][ind(0,0,0)].m[0], Hexch(0, 0, data[0].get_nb(0, 7), 0));
	}
	ftvals.close(); ftvals = File("%tvals.dat", "w", path);
	ftvals("#:t M Mx My Mz M2x M2y M2z W Wexch Wext Wanis Q1 Q2 Q3 Q4 eta eta2 eta3 eta4 PHIx PHIy PHIz THETAx THETAy THETAz  XIxx XIyy XIzz XIxy XIxz XIyz"
		   "  eta_k2 eta2_k2 eta3_k2 eta4_k2  eta_k3 eta2_k3 eta3_k3 eta4_k3   eta_k4 eta2_k4 eta3_k4 eta4_k4  Psi U_CMD UM_LL zeta  dot_eta T_sc Hz T\n");
	if(corr_max){
		corr_fout.close(); corr_fout = File("%corr.dat", "w", path);
		corr_fout.printf("#:t  eta1 eta1_2 eta1_3 eta1_4"); for(int i=0; i<corr_max; i++) corr_fout("    eta% eta%_2 eta%_3 eta%_4", i+2, i+2, i+2, i+2);
		corr_fout.printf("\n");
	}
	calc_av();	drop_tvals();
}
void Model::drop_tvals(){
	double mm = M*M, zeta = mm<1? (eta[0]-M*M)/(1-M*M): 0;
	ftvals("% %        %  %   %  %  %    %    %      %   %       %       %       %    %   %  %   %  %  % %\n",
		   t, M.abs(), M, M2, W, Q, eta, PHI, THETA, XI, eta_k2, eta_k3, eta_k4, Psi, Upsilon(M.abs(), eta[0]), UpsilonM.abs(), zeta, dot_eta, T_sc, Hext[2], T).flush();
	if(corr_max){
		corr_fout("%  %", t, eta);  for(int i=0; i<corr_max; i++) corr_fout("    %", corr[i]);
		corr_fout.printf("\n");  corr_fout.flush();
	}
}
//------------------------------------------------------------------------------
void Model::calc(int steps){
	if(!init_conditions) start_gauss();  // флаг задания н.у.
#ifndef MAGNONS
	// float sghT = sqrt(2*dt*alpha*std::max(0., std::min(T, 2*T-(T_sc+T_sc_old)/2))); // v3 2*sqrt(dt*alpha*T);
	// float sghT = sqrt(2*dt*alpha*std::max(0., std::min(T, 2*T-T_sc))); // v2 2*sqrt(dt*alpha*T);
	// float sghT = sqrt(2*dt*alpha*std::max(0., 2*T-T_sc)); // v1 2*sqrt(dt*alpha*T);
	// float sghT = sqrt(2*dt*alpha*(patchT && !calc_eq? T_sc: T)); // v4,5 2*sqrt(dt*alpha*T);
	float sghT = sqrt(2*dt*alpha*(patchT ? T_sc: T)); // v4,5 2*sqrt(dt*alpha*T);
#endif // MAGNONS
	size_t data_sz = data[0].size();  eta_old = eta[0];  float Kx2 = 2*K;
	for(int Nt=0; Nt<steps; Nt++){
		double t0 = omp_get_wtime();
		eta_old = 0;
#pragma omp parallel for reduction(+:eta_old) if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(0, i, nb, k);
				eta_old += Hex*m0;
				// WOUT(i, k, Hex); 
				Vecf<3> H = Hex + Hext + nK*m0*Kx2*nK; 
				Vecf<3> dmdt = m0%(-gamma*H -alpha*m0%H);
				m1 = m0 + .5f*dt*dmdt;
				dm = .5f*dmdt;
			}
			// exit(0);
		}
		eta_old /= data_sz*nb_sz*cell_sz;
		//-----------
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &m2 = data[2][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(1, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m1*Kx2*nK; 
				Vecf<3> dmdt = m1%(-gamma*H -alpha*m1%H);
				m2 = m0 + .5f*dt*dmdt;
				dm += dmdt;
			}
		}
		//------------
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &m2 = data[2][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(2, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m2*Kx2*nK; 
				Vecf<3> dmdt = m2%(-gamma*H -alpha*m2%H);
				m1 = m0 + dt*dmdt;
				dm += dmdt;
			}
		}
		//-------------
#ifdef MAGNONS // MAGNONS
		// делаем магнон
		/*
		Vecf<3> km(random()%nm_max, random()%nm_max, random()%nm_max); km *= 2*M_PI/(1<<data_rank);
		float phi0 = rand_alpha*random()*2*M_PI;
		Vecf<3> nm = sph_cell(rand()%sph_cells_num(5), 5), nm_perp = perp(nm);
		*/
		// MagnonsStochSrcV3 mg(data_rank, mg_split, dt*T*alpha);
		MagnonsStochField mg = MG();
#endif // MAGNONS
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			int thID = omp_get_thread_num();
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(1, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m1*Kx2*nK; 
				Vecf<3> dmdt = m1%(-gamma*H -alpha*m1%H);
				m0 += dt*_3*(dm + .5f*dmdt);
				// m0 = gauss_rotate(m0/m0.abs(), sghT, seed); // old variant, errro
				m0 /= m0.abs(); // data[2][i].m[k] = m0;

				// Vecf<3> Hm = rotate(nm_perp, nm, sin(phi0+km*coord[k]+km*zoff2pos<3>(i, data_rank)));
				
#ifndef MAGNONS
				m0 = rotate(m0, randN01.V<3>(thID)*sghT);  // ??? поворот или приращение ???
#else  // MAGNONS
				// m0 = rotate(m0, Hm*sghT);
				// m0 = rotate(m0, mg(zoff2pos<3>(i, data_rank)+coord[k])*zeta+(1-zeta)*rand_gaussV<3, float>(seed)*sghT);
				m0 = rotate(m0, mg(zoff2pos<3>(i, data_rank)+coord[k]));
#endif // MAGNONS
			}
		}
		//---------------
		/*
		unsigned int seed = 0;	
#pragma omp parallel for firstprivate(seed) if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
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
		rt_calc += omp_get_wtime() - t0;
		t += dt; if(out_tm0) ftm("% % %\n", t, data[0][0].m[0], Hexch(0, 0, data[0].get_nb(0, 7), 0));
		if(calc_eq || Nt==steps-1) calc_av();
	}
	if(calc_eq){
		fftw.fwd([&](const Ind<3> &pos, int c){ return data[0][pos].m[0][c]-M[c]; });
		std::vector<float> sp; spectr_f_max = 0; 
		fftw.spectrum(sp, spectr_f_max, 7);
		if(!spectr_eq.size()) spectr_eq.resize(sp.size(), 0.);
		for(int i=0, ssz=sp.size(); i<ssz; i++) spectr_eq[i] += sp[i];
		spectr_eq_count++;
	}
	drop_tvals();
}
//------------------------------------------------------------------------------
void Model::calc_quant(int steps){
	if(!init_conditions) start_gauss();  // флаг задания н.у.
	float sghT = sqrt(2*dt*alpha*T); // v4,5 2*sqrt(dt*alpha*T);
	size_t data_sz = data[0].size();  eta_old = eta[0];  float Kx2 = 2*K;
	for(int Nt=0; Nt<steps; Nt++){
		double t0 = omp_get_wtime();
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(0, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m0*Kx2*nK; 
				Vecf<3> dmdt = -gamma*m0%H;
				m1 = m0 + .5f*dt*dmdt;
				dm = .5f*dmdt;
			}
		}
		//-----------
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &m2 = data[2][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(1, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m1*Kx2*nK; 
				Vecf<3> dmdt = -gamma*m1%H;
				m2 = m0 + .5f*dt*dmdt;
				dm += dmdt;
			}
		}
		//------------
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &m2 = data[2][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(2, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m2*Kx2*nK; 
				Vecf<3> dmdt = -gamma*m2%H;
				m1 = m0 + dt*dmdt;
				dm += dmdt;
			}
		}
		//-------------
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k], &dm = data[3][i].m[k];
				Vecf<3> Hex = Hexch(1, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m1*Kx2*nK; 
				Vecf<3> dmdt = -gamma*m1%H;
				m0 += dt*_3*(dm + .5f*dmdt);
			}
		}
		//-------------
		float ds = 2./(1-n_s), _ds = 1./ds, dsT = ds/T;
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			int thID = omp_get_thread_num();
			ZCubeNb<3> nb = data[0].get_nb(i, 7);
			for(int k=0; k<cell_sz; k++){
				Vecf<3> &m0 = data[0][i].m[k], &m1 = data[1][i].m[k];
				Vecf<3> Hex = Hexch(0, i, nb, k);
				Vecf<3> H = Hex + Hext + nK*m0*Kx2*nK;  
				m1 = m0 - gamma*alpha*dt*(m0%(m0%H));  
				m1 = rotate(m1/m1.abs(), randN01.V<3>(thID)*sghT);
				float H_abs = H.abs();
				if(H_abs>1e-6){
					Vecf<3> nH = H/H_abs;
					/*   v0
					float p0 = m1*nH, p = roundf((1+p0)*_ds)*ds-1;
					*/
					float p0 = m1*nH, p1 = floorf((1+p0)*_ds)*ds-1, p2 = p1+ds, p;
					/*   v1
					if(p1<-1-.5f*ds) p = p2;
					else if(p2>1+.5f*ds) p = p1;
					else p = randN01.u(thID)<(p2-p0)*_ds? p1: p2; 
					*/
					if(p1<-1-.5f*ds){ p1 = -1; p2 = ds-1; }
					else if(p2>1+.5f*ds){ p1 = 1-ds; p2 = 1; }
					p = randN01.u(thID)*(1+exp(H_abs*dsT))<1? p1: p2;  // v2
					if(fabsf(p)<.999){
						//  m1 = (nH*p) + w*m_ort ==>  (nH*p)^2 + w^2*m_ort^2=1 ==> w = sqrt((1-(nH*p)^2)/m_ort^2)
						Vecf<3> m_pr = nH*p, m_ort = m1-nH*p0;
						float w = sqrtf((1-m_pr*m_pr)/(m_ort*m_ort));
						m1 = m_pr + w*m_ort; 			
					} else m1 = p*nH;
				}				
			}
		}
		for(size_t i=0; i<data_sz; ++i)	for(int k=0; k<cell_sz; k++) data[0][i].m[k] = data[1][i].m[k];
		//---------------
		rt_calc += omp_get_wtime() - t0;
		t += dt; if(out_tm0) ftm("% % %\n", t, data[0][0].m[0], Hexch(0, 0, data[0].get_nb(0, 7), 0));
		if(calc_eq || Nt==steps-1) calc_av();
	}
	if(calc_eq){
		fftw.fwd([&](const Ind<3> &pos, int c){ return data[0][pos].m[0][c]-M[c]; });
		std::vector<float> sp; spectr_f_max = 0; 
		fftw.spectrum(sp, spectr_f_max, 7);
		if(!spectr_eq.size()) spectr_eq.resize(sp.size(), 0.);
		for(int i=0, ssz=sp.size(); i<ssz; i++) spectr_eq[i] += sp[i];
		spectr_eq_count++;
	}
	drop_tvals();
}
//------------------------------------------------------------------------------
void Model::calc_av(){  // считаем средние значения
	double t0 = omp_get_wtime();
	size_t data_sz = data[0].size(); 
	double Mx = 0, My = 0, Mz = 0, M2x = 0, M2y = 0, M2z = 0, We = 0, Wa = 0, eta1 = 0, eta2 = 0, eta3 = 0, eta4 = 0,
		phi_x = 0, phi_y = 0, phi_z = 0, th_x = 0, th_y = 0, th_z = 0,
		xi_xx = 0, xi_yy = 0, xi_zz = 0, xi_xy = 0, xi_xz = 0, xi_yz = 0, psi = 0, UMx = 0, UMy = 0, UMz = 0;

	for(int &v: f_buf) v = 0;
	for(int &v: fz_buf) v = 0;
	for(float &v: f_eta_buf) v = 0;
	float _dz = fz_sz/2., _deta = f_eta_sz/2., df_eta = 1./(2./f_eta_sz*data_sz*cell_sz*nb_sz/2); 

	for(double &x: Q_buf) x = 0; 
	for(double &x: eta_k_buf) x = 0;
	for(Vec<4> &x: corr_buf) x = vec(0.);

	bool calc_Ms = Ms_start>=0 && Ms_start<data_rank;
	
#pragma omp parallel for reduction(+:Mx,My,Mz,M2x,M2y,M2z,We,Wa,eta1,eta2,eta3,eta4,phi_x,phi_y,phi_z,th_x,th_y,th_z,xi_xx,xi_yy,xi_zz,xi_xy,xi_xz,xi_yz,psi,UMx,UMy,UMz) if(threads!=1)
	for(size_t cID=0; cID<data_sz; cID++){  // начало цикла по ячейкам
		int thID = omp_get_thread_num();
		ZCubeNb<3> nb = data[0].get_nb(cID, 7);
		Vecf<3> &ms = Ms_arr[0][data[0].offset2pos(cID)]; if(calc_Ms) ms = Vecf<1>();
		for(int k=0; k<cell_sz; k++){  // начало цикла внутри ячейки
			const Vecf<3> &m0 = data[0][cID].m[k]; if(calc_Ms) ms += m0;
			Mx += m0[0];        My += m0[1];        Mz += m0[2];
			M2x += m0[0]*m0[0]; M2y += m0[1]*m0[1]; M2z += m0[2]*m0[2];
			We -= Hext*m0;
			float nKm = nK*m0; Wa -= K*nKm*nKm;
			Vecf<3> MnK = m0%nK, phi = MnK*nKm, th = m0%MnK*nKm;
			phi_x += phi[0]; phi_y += phi[1]; phi_z += phi[2];
			th_x  -=  th[0];  th_y -=  th[1];  th_z -=  th[2];

			xi_xx -= m0[0]*m0[0]-1;  
			xi_yy -= m0[1]*m0[1]-1;  
			xi_zz -= m0[2]*m0[2]-1;  
			xi_xy -= m0[0]*m0[1]; 
			xi_xz -= m0[0]*m0[2]; 
			xi_yz -= m0[1]*m0[2]; 

			if(f_rank>=0) f_buf[thID*f.size()+f.find(m0)]++;
			if(fz_sz){ int fid = floor((1+m0[2])*_dz); fz_buf[thID*fz_sz+(fid<0 ? 0 : (fid>=fz_sz ? fz_sz-1 : fid))]++; }

			for(int l=0; l<corr_max; l++) for(int d=0; d<corr_direct_sz; d++){  // расчет корреляционной функции
					float e = m0*data[0][corr_direct_offs[(cID*corr_max + l)*corr_direct_sz + d]].m[k];
					Vec<4> &dste = corr_buf[thID*corr_max+l];  dste[0] += e; float ee = e; for(int i=1; i<4; i++){ ee *= e; dste[i] += ee; }
				} 
#ifdef CALC_Q
			Vecf<3> mj_buf[nb_sz]; int j = 0;  // массив магнитных моментов соседей
#endif  // CALC_Q
			for(const auto &l: nb_pos[k]){  // начало цикла по соседям
				const Vecf<3> &mj = data[0][cID+nb[l.off]].m[l.cell];
				float e = m0*mj, e2 = e*e; eta1 += e; eta2 += e2;  eta3 += e*e2; eta4 += e2*e2; 
				psi -= m0*(mj%(mj%nK))*(mj*nK);
				Vecf<3> UM = m0%(m0%mj); UMx += UM[0]; UMy += UM[1]; UMz += UM[2];

				if(f_eta_sz){ int fid = floor((1+e)*_deta); f_eta_buf[thID*f_eta_sz+(fid<0 ? 0 : (fid>=f_eta_sz ? f_eta_sz-1 : fid))] += df_eta; }
#ifdef CALC_Q
				mj_buf[j++] = mj;  // накопление соседей
#endif  // CALC_Q
			}  // конец цикла по соседям

			//------- собственно расчет Q ------------
#ifdef CALC_Q  
			int qt = 0, ssz = Qtable_sz[k][0];  j = 0;  // общий счетчик в Qtable
			for(const Ind<2>& ik: Qtable[k]){
				float e = mj_buf[ik[0]]*mj_buf[ik[1]], en = e, Qk = mj_buf[ik[0]]*(m0%(m0%mj_buf[ik[1]]));
				int off = (thID*Q_sz+qt)*4; eta_k_buf[off] += e; for(int i=1; i<4; i++){ en *= e; eta_k_buf[off+i] += en; }
				Q_buf[thID*Q_sz+qt] += Qk;
				if(++j>=ssz){ qt++; if(qt<Q_sz) ssz += Qtable_sz[k][qt]; }
			}
#endif  // CALC_Q
			//------- конец расчета Q ----------------- 
		}  // конец цикла внутри ячейки
	} // конец цикла по ячейкам

	M  = vec(Mx, My, Mz)/(data_sz*cell_sz);  
	M2 = vec(M2x, M2y, M2z)/(data_sz*cell_sz);
	W  = vec(-J*eta1/2+We+Wa, -J*eta1/2, We, Wa)/(data_sz*cell_sz);
	eta = vecf(eta1, eta2, eta3, eta4)/(data_sz*cell_sz*nb_sz);

	PHI = vec(phi_x, phi_y, phi_z)/(data_sz*cell_sz);
	THETA = vec(th_x, th_y, th_z)/(data_sz*cell_sz);
	XI = vec(xi_xx, xi_yy, xi_zz, xi_xy, xi_xz, xi_yz)/(data_sz*cell_sz);
	
	Psi = psi/(data_sz*cell_sz*nb_sz);
	UpsilonM = -vec(UMx, UMy, UMz)/(2*data_sz*cell_sz*nb_sz); 

	for(Vec<4> &x: corr) x = vec(0.);
	for(int i=0, sz=corr_buf.size(); i<sz; i++) corr[i%corr_max] += corr_buf[i]/(data_sz*cell_sz*corr_direct_sz);
	
#ifdef CALC_Q
	Q = eta_k2 = eta_k3 = eta_k4 = vec(0.);
	for(int i=0, sz=Q_buf.size(); i<sz; i++) Q[i%Q_sz] += Q_buf[i]/(data_sz*cell_sz*Qtable_sz[0][i%Q_sz]);
	for(int i=0, sz=eta_k_buf.size(); i<sz; i++){
#ifdef FCC
		int qt = i/4%Q_sz-1, qtt = qt+1; if(qt<0) continue;
#else  // FCC
		int qt = i/4%Q_sz, qtt = qt;
#endif // FCC
		(qt==0? eta_k2: (qt==1? eta_k3: eta_k4))[i%4] += eta_k_buf[i]/(data_sz*cell_sz*Qtable_sz[0][qtt]);
	}
	// Z2 z2; Z2.calc(M.abs(), eta);
#endif // CALC_Q

	if(calc_Ms)	for(int R=1; R<data_rank; R++){
			Mesh<Vecf<3>, 3> &src = Ms_arr[R-1], &dst = Ms_arr[R]; dst.fill(Vecf<3>()); 
			for(const Ind<3>& pos: irange(src.bbox())) dst[pos/2] += src[pos];
		}
	dot_eta = (eta[0]-eta_old)/dt;
#ifdef CALC_Q
	// T_sc_old = T_sc;
	// if(patchT) T_sc = 2*(-dot_eta/(4*alpha) + Hext*UpsilonM + K*Psi  - .5*J*(eta[1]-1+Qcoeff*Q))/(eta[0]+eta_old); 	
	T_sc =  (1-T_sc_weight)*T_sc + T_sc_weight*std::max(0., 2*T - (-dot_eta/(4*alpha) + Hext*UpsilonM + K*Psi  - .5*J*(eta[1]-1+Qcoeff*Q))/eta[0]);
#endif //CALC_Q
	if(calc_eq){ 
		Meq += M; Mabs_eq += M.abs(); M2eq += M2; Weq += W; Qeq += Q; eta_eq += eta;
		PHIeq += PHI; THETAeq += THETA; XIeq += XI; Psi_eq += Psi; UpsilonMeq += UpsilonM;
		eta_k2_eq += eta_k2;  eta_k3_eq += eta_k3;  eta_k4_eq += eta_k4;
		for(int i=0; i<corr_max; i++) corr_eq[i] += corr[i];
		
		Tsc_eq += T_sc;  MxMeq += M&M;
		// Seq += S;
	}

	if(f_rank>=0){
		float  df = 1/(4*M_PI/f.size()*data_sz*cell_sz);
		for(int i=0, sz=f_buf.size(), fsz=f.size(); i<sz; i++) f[i%fsz] += df*f_buf[i];
		if(calc_eq) for(int i=0, sz=f.size(); i<sz; i++) f_eq[i] += f[i];
	}
	if(fz_sz>0){
		for(float &v: fz) v = 0;
		float  df = 1./(2./fz_sz*data_sz*cell_sz);
		for(int i=0, sz=fz_buf.size(); i<sz; i++) fz[i%fz_sz] += df*fz_buf[i];
		if(calc_eq) for(int i=0; i<fz_sz; i++) fz_eq[i] += fz[i];
	}
	if(f_eta_sz>0){
		for(float &v: f_eta) v = 0;
		for(int i=0, sz=f_eta_buf.size(); i<sz; i++) f_eta[i%f_eta_sz] += f_eta_buf[i];
		if(calc_eq) for(int i=0; i<f_eta_sz; i++) f_eta_eq[i] += f_eta[i];
	}
	if(calc_eq) eq_count++;
	rt_diagn += omp_get_wtime() - t0;
}
//------------------------------------------------------------------------------
void Model::clean_av_eq(){
	eq_count = 0; Meq = vec(0.); Mabs_eq = 0.; M2eq = vec(0.); Weq = vec(0.); Qeq = vec(0.);
	eta_eq = vec(0.); eta_k2_eq = vec(0.); eta_k3_eq = vec(0.); eta_k4_eq = vec(0.);
	PHIeq = vec(0.); THETAeq = vec(0.); XIeq = vec(0.); Psi_eq = 0; UpsilonMeq = Vec<3>();
	MxMeq = vec(0.);
	// UpsHextMeq = 0;   HextMMMeq = 0;

	if(f_rank>=0) f_eq.fill(0.f);
	for(float &v: fz_eq) v = 0;
	for(float &v: f_eta_eq) v = 0;
	Tsc_eq = 0; // Seq = 0;

	spectr_eq_count = 0; spectr_eq.clear(); spectr_f_max = 0;
}
//------------------------------------------------------------------------------
void Model::finish(){
	if(eq_count){
		Weq /= eq_count;
		Meq /= eq_count;
		M2eq /= eq_count;
		Mabs_eq /= eq_count;
		Qeq /= eq_count;
		eta_eq /= eq_count;
		eta_k2_eq /= eq_count;
		eta_k3_eq /= eq_count;
		eta_k4_eq /= eq_count;
		PHIeq /= eq_count;
		THETAeq /= eq_count;
		XIeq /= eq_count;
		Psi_eq /= eq_count;
		UpsilonMeq /= eq_count;
		Tsc_eq /= eq_count;
		// Seq /= eq_count;
		
		if(f_rank>=0) for(int i=0, sz=f.size(); i<sz; i++) f_eq[i] /= eq_count;
		if(fz_sz>=0)  for(int i=0; i<fz_sz; i++) fz_eq[i] /= eq_count;
		if(f_eta_sz>=0)  for(int i=0; i<f_eta_sz; i++) f_eta_eq[i] /= eq_count;
		for(int i=0; i<corr_max; i++) corr_eq[i] /= eq_count;

		MxMeq /= eq_count;
	} else {
		Weq = W;
		Meq = M;
		M2eq = M2;
		Mabs_eq = M.abs();
		if(f_rank>=0) f_eq = f.copy();
		Qeq = Q;
		eta_eq = eta;
		eta_k2_eq = eta_k2;
		eta_k3_eq = eta_k3;
		eta_k4_eq = eta_k4;
		PHIeq = PHI;
		THETAeq = THETA;
		XIeq = XI;
		Psi_eq = Psi;
		UpsilonMeq = UpsilonM;
		Tsc_eq = T_sc;
		// Seq = S;

		if(f_rank>=0) for(int i=0, sz=f.size(); i<sz; i++) f_eq[i] = f[i];
		if(fz_sz>=0)  for(int i=0; i<fz_sz; i++) fz_eq[i] = fz[i];
		if(f_eta_sz>=0)  for(int i=0; i<f_eta_sz; i++) f_eta_eq[i] = f_eta[i];

		corr_eq = corr;

		MxMeq = M&M;
	}
	if(corr_max){
		File  fcorr_eq("%corr_eq.dat", "w", path); fcorr_eq("#:k eta eta2 eta3 eta4\n1 %\n", eta_eq);
		for(int k=0; k<corr_max; k++) fcorr_eq("% %\n", k+2, corr_eq[k]);		
	}
	//-------------------------
	if(spectr_eq_count) for(auto& s: spectr_eq) s /= spectr_eq_count;
	else{
		fftw.fwd([&](const Ind<3> &pos, int c){ return data[0][pos].m[0][c]-M[c]; });
		std::vector<float> sp; spectr_f_max = 0; 
		fftw.spectrum(sp, spectr_f_max, 7);
		if(!spectr_eq.size()) spectr_eq.resize(sp.size(), 0.);
		for(int i=0, ssz=sp.size(); i<ssz; i++) spectr_eq[i] = sp[i];	
	}
	std::ofstream fout(path+"spectr-eq.dat");  fout<<"#:nu ampl\n";
	for(int i=0, ssz=spectr_eq.size(); i<ssz; i++) fout<<(i+.5)*spectr_f_max/ssz<<' '<<spectr_eq[i]<<'\n';	
}
//------------------------------------------------------------------------------
void Model::dump_fz(const char *path, bool eq) const {
	std::ofstream fout(path); fout<<"#:m_z f\n";
	for(int i=0; i<fz_sz; i++) fout<<(i+.5)*2/fz_sz-1<<' '<<(eq? fz_eq: fz)[i]<<'\n';
}
void Model::dump_f_eta(const char *path, bool eq) const {
	std::ofstream fout(path); fout<<"#:eta f\n";
	for(int i=0; i<f_eta_sz; i++) fout<<(i+.5)*2/f_eta_sz-1<<' '<<(eq? f_eta_eq: f_eta)[i]<<'\n';
}
void Model::dump_Ms_arr(){
	if(Ms_start<0 || Ms_start>=data_rank) return;
	if(!fMs.size()){ fMs.resize(Ms_arr.size()-Ms_start); for(int i=Ms_start, sz=Ms_arr.size(); i<sz; i++) fMs[i].open("%/Ms%.msh", "w", path, i); }
	for(int i=Ms_start, sz=Ms_arr.size(); i<sz; i++) Ms_arr[i].dump(fMs[i]);
}
//------------------------------------------------------------------------------
void Model::dump_data(const char *path) const {
	char head[64]; for(int i=0; i<64; i++){ head[i] = 0; } snprintf(head, 64, "%s %i", mode, data_rank);
	File fout(path, "w"); fout.write(head, 64); fout.write(&(data[0][0]), data[0].size()*sizeof(Cell));
}
bool Model::load_data(const char *path){
	double t0 = omp_get_wtime();
	File fin(path, "r"); char head[64];  if(fin.read(head, 64)!=64) return false;
	std::string rmode; int rrank; std::stringstream(head)>>rmode>>rrank;
	if(rmode!=mode){ WMSG(mode, rmode, data_rank, rrank);  return false; }
	data_rank = rrank; for(int k=0; k<4; k++) data[k].init(data_rank);
	Natoms = data[0].size()*cell_sz;
	if(fin.read(&(data[0][0]), data[0].size()*sizeof(Cell))!=data[0].size()*sizeof(Cell)){ WMSG("file too short"); return false; }
	rt_init += omp_get_wtime()-t0;
	open_tvals(); init_conditions = true;  // флаг задания н.у.
	return true;
}
//------------------------------------------------------------------------------
void Model::check_rand(int steps, const char *path){  ///< создает файл с приращениями <M>, <eta> от случайного источника по числу шагов, остальные члены отключены.
	float sghT = sqrt(2*dt*alpha*T); 
	std::ofstream fout(path);  fout<<"#:dM dMx dMy dMz deta\n";
	size_t data_sz = data[0].size(); 
	for(int it=0; it<steps; it++){
#pragma omp parallel for if(threads!=1)
		for(size_t i=0; i<data_sz; ++i){
			int thID = omp_get_thread_num();
			for(int k=0; k<cell_sz; k++){
				const Vecf<3> &m0 = data[0][i].m[k]; Vecf<3> &m1 = data[1][i].m[k];
				m1 = rotate(m0, randN01.V<3>(thID)*sghT);
			}
		}
		double M1x=0, M1y=0, M1z=0, eta1=0;
	
#pragma omp parallel for reduction(+:M1x,M1y,M1z,eta1) if(threads!=1)
		for(size_t cID=0; cID<data_sz; cID++){  // начало цикла по ячейкам
			// int thID = omp_get_thread_num();
			ZCubeNb<3> nb = data[0].get_nb(cID, 7);
			for(int k=0; k<cell_sz; k++){  // начало цикла внутри ячейки
				const Vecf<3> &m1 = data[1][cID].m[k];
				M1x += m1[0];        M1y += m1[1];        M1z += m1[2];
				for(const auto &l: nb_pos[k]){  // начало цикла по соседям
					const Vecf<3> &mj = data[1][cID+nb[l.off]].m[l.cell];
					float e = m1*mj; eta1 += e; 
				}  // конец цикла по соседям
				//------- конец расчета Q ----------------- 
			}  // конец цикла внутри ячейки
		} // конец цикла по ячейкам

		Vec<3> M1 = vec(M1x, M1y, M1z)/(data_sz*cell_sz);  
		eta1 /= (data_sz*cell_sz*nb_sz);

		fout<<M1.abs()-M.abs()<<' '<<M1-M<<' '<<eta1-eta[0]<<'\n';		
	}
}
//------------------------------------------------------------------------------
void Model::calc_spectrum(const char *path){
	fftw.fwd([&](const Ind<3> &pos, int c){ return data[0][pos].m[0][c]-M[c]; });
	std::vector<float> sp; float f_max = 0; 
	fftw.spectrum(sp, f_max, 7);
	if(path){
		std::ofstream fout(path);  fout<<"#:nu ampl\n";
		for(int i=0, sz=sp.size(); i<sz; i++) fout<<(i+.5)*f_max/sz<<' '<<sp[i]<<'\n';
	} else {		
		for(int i=0, sz=sp.size(); i<sz; i++) fspectrum<<t<<' '<<(i+.5)*f_max/sz<<' '<<sp[i]<<'\n';
		fspectrum<<std::endl;
	}
}
//------------------------------------------------------------------------------
