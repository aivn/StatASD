#include "magnons.hpp"
//------------------------------------------------------------------------------
MagnonsStochSrc::MagnonsStochSrc(int data_rank, int mg_split, float hTalpha){
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
//------------------------------------------------------------------------------
MagnonsStochSrcV2::MagnonsStochSrcV2(int data_rank, int mg_split, float hTalpha){
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
//------------------------------------------------------------------------------
MagnonsStochSrcV0::MagnonsStochSrcV0(int data_rank, int mg_split, float hTalpha){
	for(int i=0; i<3; i++) k[i] = 2*M_PI/(1<<data_rank)*(rand()%(mg_split+1));
	A = sqrt(2*hTalpha)*rand_gauss();
	sph_init_table(5);
	n = sph_cell(rand()%sph_cells_num(5), 5);
	n_perp = perp(n);
	phi0 = rand()*rand_alpha*2*M_PI;
}
//------------------------------------------------------------------------------
MagnonsStochSrcV3::MagnonsStochSrcV3(int data_rank, int mg_split, float hTalpha){
	sph_init_table(5);
	for(int i=0; i<3; i++) k[i] = 2*M_PI/(1<<data_rank)*(rand()%(mg_split+1));
	A = sqrt(2*hTalpha)*rand_gauss();
	n = sph_cell(rand()%sph_cells_num(5), 5);
	n_perp = perp(n);
	phi0 = rand()*rand_alpha*2*M_PI;
	
	magnons.init(Ind<3>(mg_split), Vecf<3>(), vecf(1<<data_rank));
	magnons.periodic = 7;
	magnons.interp = 0x333;
	for(auto &mg: magnons) mg = 4*sph_cell(rand()%sph_cells_num(5), 5);
}
//------------------------------------------------------------------------------
MagnonsStochField::MagnonsStochField(int data_rank, Ind<3> mg_sz, float A_){
	sph_init_table(5); A = A_;
	k = mg_sz*(2*M_PI/(1<<data_rank)); 
	n = sph_cell(rand()%sph_cells_num(5), 5);
	n_perp = perp(n);
	phi0 = rand()*rand_alpha*2*M_PI;
	/*	
	magnons.init(mg_sz, Vecf<3>(), vecf(1<<data_rank));
	magnons.periodic = 7;
	magnons.interp = 0x333;
	for(auto &mg: magnons) mg = 4*sph_cell(rand()%sph_cells_num(5), 5);
	*/
}
//------------------------------------------------------------------------------
const int n_b = 8;
const Vecf<3> stencil[4] = {
	vecf(-.5f, -.5f,  .5f),
	vecf(-.5f,  .5f,  .5f),
	vecf( .5f, -.5f,  .5f),
	vecf( .5f,  .5f,  .5f)
};
const float hbar = .5f;

void MagnonsStochFieldGen::init(int data_rank_, float T, float cL, float halpha_){
	data_rank = data_rank_; halpha = halpha_;
	int L = 1<<data_rank; float _N = 1./(2*L*L*L);
	// for(Ind<3> pos: irange<3>(L)){
	// for(Ind<3> pos: irange<3>(-L/2-1, L/2+2)){ // MG14
	for(Ind<3> pos0: irange<3>(-L, L+1)){ // MG15
		Ind<3> pos(abs(pos0[0]), abs(pos0[1]), abs(pos0[2]));
		// if(pos.abs()<.1*L) continue; // MG15
		// if(pos.abs()>.61*L) continue; // MG16
		if(!pos[0] && !pos[1] && !pos[2]) continue;
		Vecf<3> k = 2*M_PI/L*pos;
		bool k_ok = true; float sigma = 0;
		for(const Vecf<3> &d: stencil){
			float dphi = k*d; if(dphi>M_PI){ k_ok = false; break; }
			sigma += 1-cosf(dphi);
		}
		if(!k_ok) continue;
		// float ek = hbar*sigma, fk = 24.f/(cL*(expf(ek/T)-1));
		// float ek = hbar*(.5+sigma), fk = 24.f/(cL*(expf(ek/T)-1)); // MG13
		// float ek = hbar*(sigma), fk = 1./(expf(ek/T)-1); // MG14
		// float ek = hbar*(sigma), fk = 1./(expf(ek/T)-1); // MG15
		// float ek = hbar*(sigma), fk = .25/(expf(ek/T)-1); // MG17
		float ek = hbar*(sigma), fk = 1./(expf((.25+ek)/T)-1); // MG18
		node_t node; node.k = pos; 
		// node.A = sqrtf(hbar*fk*_N);
		node.A = 8*sqrtf(fk*ek*_N)/M_PI;  // MG0, MG6 (см. хидер), MG13, MG14, MG15
		// node.A = 8*sqrtf(fk*hbar*_N)/M_PI;   // MG1
		// node.A = 8*sqrtf(fk*ek*sigma*_N)/M_PI;   // MG2
		// MG3 и далее основаны на неравных вероятностях k
		// node.P = sigma; node.A = 8*sqrtf(fk*ek*_N)/M_PI;  // MG3
		// node.P = sigma; node.A = 8*sqrtf(fk*hbar*_N)/M_PI;  // MG4
		// node.P = sqrt(sigma); node.A = 8*sqrtf(fk*ek*_N)/M_PI;  // MG5
		// float ssig = sqrt(sigma); node.A = 8*sqrtf(fk*ek*_N/ssig)/M_PI;  // MG7
		// float ssig = pow(sigma, 1./3); node.A = 8*sqrtf(fk*ek*_N/ssig)/M_PI;  // MG8
		// float ssig = pow(sigma, .25); node.A = 8*sqrtf(fk*ek*_N/ssig)/M_PI;  // MG9
		// float ssig = pow(sigma, .5); node.A = 8*sqrtf(fk*ek*_N*ssig)/M_PI;  // MG10
		// float ssig = pow(sigma, .25); node.A = 8*sqrtf(fk*ek*_N*ssig)/M_PI;  // MG11
		// float ssig = pow(sigma, .125); node.A = 8*sqrtf(fk*ek*_N*ssig)/M_PI;  // MG12
		
		table.push_back(node);
	}
	WMSG(table.size());
	float sqrtNk = sqrtf(table.size());
	Pmax = 0; for(auto &n: table){ Pmax += n.P; n.P = Pmax; n.A *= sqrtNk; }
}
//------------------------------------------------------------------------------
/*  MG0...2 */
float MagnonsStochFieldGen::get_A(Ind<3> &k) const {
	const node_t &node = table[int(rand_alpha*rand()*table.size())];
	k = node.k; return node.A*sqrtf(halpha);
}
/*
float MagnonsStochFieldGen::get_A(Ind<3> &k) const {
	float P = Pmax*rand_alpha*rand();
	int ai = 0, bi = table.size();
	while(bi-ai>1){  // дихотомия, bi не включается, ответ будет в ai
		int ci = (ai+bi)/2;
		if(table[ci].P<P) ai = ci;
		else if(ci==0 || table[ci-1].P<P){ ai = ci; break; } // точное попадание
		else bi = ci;
	}
	k = table[ai].k; return table[ai].A;
	}*/
MagnonsStochField MagnonsStochFieldGen::operator () () const {
	Ind<3> k; float A = get_A(k);
	return MagnonsStochField(data_rank, k+ind(1), A);
}
//------------------------------------------------------------------------------
