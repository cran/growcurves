#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmCmiar(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Hmat, SEXP Wcase, SEXP Wperson,
        SEXP Omega, SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec,
        SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd, SEXP corsessInt, SEXP typeMM)
{
BEGIN_RCPP
    // Time run
    //wall_clock timer;
    //timer.tic();
    // Wrap SEXP objects in Rcpp Vector and Matrix objects.  Not copies,
    // but wrappers.
    // Data: observations, fe, re, mm weights, car weights, session groups,
    // client id's
    // Counters (number of); cbt cases, mcmc iter, burnin.
    NumericVector yr(yvec);
    NumericMatrix Xr(Xmat);
    NumericMatrix Zr(Zmat);
    NumericMatrix Hr(Hmat);
    NumericMatrix Wcr(Wcase);
    NumericMatrix Wp(Wperson);
    NumericMatrix Or(Omega);
    NumericVector opr(omegaplusvec);
    IntegerVector gr(groupsvec);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    int corsess = as<int>(corsessInt); /* correlation between sets of sess eff*/
    int typemm = as<int>(typeMM); /* typemm == 0 for ind prior on umat, 1 for CAR prior on umat */
    double ustrength = as<double>(ustrengthd);
    int nkeep = (niter -  nburn)/nthin;

    // Extract key row and column dimensions we will need to sample parameters
    int nc = Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), nv = Hr.ncol(), ns = Or.ncol();
    int npcbt = Wp.nrow(); int nrho = 0.5*nv*(nv-1);

    // Compute np, number of unique clients
    // algorithm assumes cases clustered by unique client
    IntegerVector dpr = diff(pr);
    icolvec diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    /* int np = accumulate(dpr.begin,dpr.end,1); */

    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
    mat hmat(Hr.begin(), nc, nv, false);
    mat wcase(Wcr.begin(),nc, ns, false);
    mat wpers(Wp.begin(),npcbt, ns, false); /* mmwt */
    mat omega(Or.begin(), ns, ns, false);
    colvec y(yr.begin(), nc, false);
    colvec omegaplus(opr.begin(), ns, false);
    icolvec groups(gr.begin(), ns, false);
    icolvec persons(pr.begin(), nc, false);

    // pre-compute quadratic product of design matrix for use in sampling
    field<mat> zsplit(np,1); field<mat> ztzsplit(np,1);
    field<mat> xsplit(np,1); field<mat> wsplit(np,1);
    field<mat> hsplit(np,1);
    field<colvec> ysplit(np,1); 
    CovUnitSplitMV(zmat, xmat, wcase, hmat, y, zsplit, ztzsplit, xsplit, wsplit, hsplit, ysplit, persons);
    mat xtx; prodMatOne(xmat,xtx); 
    field<mat> hmatws(ns,1), hwsthws(ns,1);
    prodHW(wcase, hmat, hmatws, hwsthws); /* multiple each nc x nv, hmat matrix by column s of wcase, ws - data object to use for sampling MM effects  */

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    //arma_rng::set_seed_random(); /* arma */

    // Armadillo structures to capture results
    // Will be wrapped into a list element on return
    colvec Deviance(nkeep); Deviance.zeros();
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    colvec Alpha(nkeep);
    mat Beta(nkeep,nf);
    mat B(nkeep,nr*np);
    rowvec brow(np*nr); brow.zeros();
    mat MM(nkeep,nv*npcbt);
    rowvec mmrow(npcbt*nv); mmrow.zeros();
    mat U(nkeep,nv*ns);
    rowvec urow(ns*nv); urow.zeros();
    mat Taub(nkeep,nr);
    mat Tauu(nkeep,nv);
    rowvec taurow(nv); taurow.zeros();
    mat Rhotauu(nkeep,nrho);
    rowvec rhotaurow(nrho); rhotaurow.zeros();
    //colvec Taubeta(nkeep);
    //colvec Taualph(nkeep);
    colvec Taue(nkeep);

    // Initialize parameter values
    double taubeta, taualph, taue;
    taubeta = taualph = taue = 1;
    rowvec taub(nr); taub.ones();
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    mat    bmat = randn<mat>(np,nr)*sqrt(1/taub(0));
    // capture session effects in ns x nsr matrix, umat
    mat S(nv,nv); /* variance of u */
    mat L(nv,nv); int i, j; /* precision for umat */
    for(i = 0; i < nv; i++)
    {
        S(i,i) = 1/ustrength;
    }
    for(i = 0; i < nv; i++)
    {
        for( j = i+1; j < nv; j++)
        {
            S(i,j) = S(j,i) = corsess*sqrt(S(i,i))*sqrt(S(j,j));
        } /*end loop j */
    } /* end loop i filling variance matrix for umat, S */
    L = inv(S); mat V(nv,nv); V = L; /* V is prior for L */
    double nu = nv; /* df for L ~ W(V,nu) */
    /* compute correlation matrix for each iteration */
    mat Linv(nv,nv); int nelem, totelem; double rho;
    mat umat(ns,nv); rowvec m(nv); m.zeros(); rmvnrnd(L,umat,m);
    mat mmmat(npcbt,nv); mmmat.zeros(); /* client effect mapped from session */
    colvec zb(nc); zb.zeros();
    int ng; /* number of groups for MM (session) effects */
    mat qmat(ns,ns);
    if( typemm == 0)  /* mmi */
    {
	qmat.eye(ns,ns);
	ng = 0;
    }else{ /* typemm = 1 */
    	qmat = -omega; qmat.diag() = omegaplus;/* prior precision for joint u */
	// compute ng, number of session groups
    	IntegerVector dgr = diff(gr);
    	icolvec diffgroups(dgr.begin(),ns - 1, false);
    	ng = accu(diffgroups) + 1;
    }
    colvec resid(nc); resid.zeros(); /* regression residual */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    /* goodness of fit related measures */
    double deviance; colvec devres(4); devres.zeros(); 
    rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0, kk;
    
    // set hyperparameter values
    double a2, a4, b2, b4;
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */
    //a5 = 0.01; b5 = 0.01; /* taualph */
    //a6 = 0.1; b6 = 0.1; /* taubeta */

    // Conduct posterior samples and store results
    int k, l, h;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) Rcout << "Interation: " << k << endl;
        bmvstep(zsplit, ztzsplit, xsplit, wsplit, hsplit, ysplit, beta, umat, alpha, taue, taub,
            zb, bmat, np, nr);
        betamvlsstep(xmat, xtx, wcase, hmat, y, beta, umat, alpha, taue, 
                    zb, nf);
	if( typemm == 0)
	{
		uindmvstep(xmat, wcase, wpers, hmat, hmatws, hwsthws, zb, beta, y, umat,
                	mmmat, alpha, taue, L, ns, nc);
	}else{ /* typemm == 1, MM-MCAR */
		umvstep(xmat, omega, wcase, wpers, hmat, hmatws, hwsthws, zb, beta, y, omegaplus, umat,
                	mmmat, alpha, taue, L, ns, nc);
	}
        alphamvlsstep(xmat, wcase, beta, hmat, zb, y, umat, resid, alpha,
                    taue, nc);
        taumvlsstep(bmat, resid, qmat, umat, V, L, taub, taue,
                nu, a2, a4, b2, b4, ns, np, nc, ng);

        if(k >= nburn)
        {
            deviance =  dev(resid, taue); /* scalar double deviance */
            dmarg(resid, taue, devmarg); /* 1 x nc vector of densities*/
            kk = k - nburn;
	    if( kk == ((oo+1)*nthin - 1) )
            {
            	Deviance(oo) = deviance;
            	Devmarg.row(oo)  = devmarg;
            	Beta.row(oo) = trans(beta);
            	Alpha(oo) = alpha;
            	for( l = 0; l < nr; l++)
            	{
                	brow.cols( np*l,(np*(l+1)-1) ) = trans( bmat.col(l) );
            	}
            	/* nv sets of session effects */
            	for( h = 0; h < nv; h++)
            	{
                	urow.cols( ns*h,(ns*(h+1)-1) ) = trans( umat.col(h) );
                	mmrow.cols( npcbt*h,(npcbt*(h+1)-1) ) = trans( mmmat.col(h) );
                	taurow.col(h) = L(h,h);
            	}
            	/* correlation coefficients between sets of session effects */
            	rho = 0; Linv = inv(L); totelem = 0; nelem = 0;
            	for( h = 0; h < (nv-1); h++)
            	{
                	nelem = nv - (h+1);
                	for(l = (h+1); l < nv; l++)
                	{
                    		rho = Linv(h,l)/( sqrt(Linv(h,h))*sqrt(Linv(l,l)) );
                    		rhotaurow.col(totelem + (l - h) - 1) = rho;
                	}
                	totelem += nelem;
            	}
            	B.row(oo) = brow;
            	U.row(oo) = urow;
            	MM.row(oo) = mmrow;
            	Tauu.row(oo) = taurow;
            	Rhotauu.row(oo) = rhotaurow;
            	//Taubeta(oo) = taubeta;
            	//Taualph(oo) = taualph;
            	Taue(oo) = taue;
            	Taub.row(oo) = taub;
            	Resid.row(oo) = trans(resid);
		oo += 1; /* increment sample counter */
            } /* end conditional statement on returning sample */
        } /* end if k > burnin, record results for return */

    } /* end MCMC loop over k */

    // DIC
    dic1mvcomp(Deviance, Alpha, Beta, B, Taue, U, xmat, zmat, hmat, wcase,
            y, persons, devres, np); /* devres = c(dic,dbar,dhat,pd) */
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  Rcpp::Named("devres") = devres,
                                  Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("U") = U,
                                  Rcpp::Named("MM") = MM,
                                  Rcpp::Named("Tauu") = Tauu,
                                  Rcpp::Named("Rhotauu") = Rhotauu,
                                  //Rcpp::Named("Taubeta") = Taubeta,
                                 // Rcpp::Named("Taualph") = Taualph,
                                  Rcpp::Named("Taue") = Taue,
                                  Rcpp::Named("Taub") = Taub,
                                  Rcpp::Named("Residuals") = Resid
				  );

    // Print CPU runtime
     //double n_secs = timer.toc();
     //cout << "took " << n_secs << " seconds" << endl;
END_RCPP
} /* end MCMC function returning SEXP */

    // static function to sample random effects
    // updates bmat (npxnr) and zb[{i:person(i)==j}] = {zj * bmat[j,]}
    SEXP bmvstep(const field<mat>& zsplit, const field<mat>& ztzsplit,
	    const field<mat>& xsplit, const field<mat>& wsplit, const field<mat>& hsplit,
            const field<colvec>& ysplit, const colvec& beta,
            const mat& umat, double alpha, double taue,
            const rowvec& taub, colvec& zb,
            mat& bmat, int np, int nr)
    {
	BEGIN_RCPP
        // Compute cholesky of prior precision, pbmat, for mvn sampling
        // of bmat through cholesky decomposition
        mat pbmat; pbmat.eye(nr,nr);
        int i, nj; int startrow = 0, endrow;
        for(i = 0; i < nr; i++)
        {
            pbmat.diag()(i) = taub(i);
        }
        mat zj, wj, hj, zjtzj, xj;
        colvec cj, yj, bj(nr);

        // sample bmat, independently, by person
        int np = zsplit.n_rows;
        int j;
        for(j=0; j < np; j++)
        {
            /* extract by-person matrix views from data*/
	    zj = zsplit(j,0); yj = ysplit(j,0); zjtzj = ztzsplit(j,0); xj = xsplit(j,0);
	    wj = wsplit(j,0); hj = hsplit(j,0); 
            /* sample posterior of nj block of bmat*/
            cj = alpha + xj*beta;
            int nv = umat.n_cols;
            for(i = 0; i < nv; i++)
            {
                cj += hj.col(i) % (wj*umat.col(i));
            }
            bj.zeros();
            rmvnchol(zj, zjtzj, pbmat, yj, cj, bj, nr, taue);
            bmat.row(j) = trans(bj);

            // compute zbj (jth piece) of zb = [z1*b1vec,...,znp*bnpvec]
	    nj 				= zsplit(j,0).n_rows;
	    endrow 			= startrow + nj - 1;
            zb.rows(startrow,endrow) 	= zj*bj;
	    startrow 			+= nj;

        } /* end loop j for sampling bj*/
        END_RCPP

    } /* end function bstep for sampling bmat and zb */

    SEXP betamvstep(const mat& xmat, const mat& xtx, const mat& wcase, const mat& hmat,
                const colvec& y, colvec& beta, const mat& umat, double alpha,
                double taue, double taubeta, const colvec& zb, int nf)
    {
        BEGIN_RCPP
        // Compute cholesky of prior precision, pbmat, for mvn sampling
        mat pbetamat; pbetamat.eye(nf,nf); pbetamat *= taubeta;

        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb; /* nc x 1 */
        int i; int nv = umat.n_cols;
        for(i = 0; i < nv; i++)
        {
            c += hmat.col(i) % (wcase*umat.col(i));
        }
        rmvnchol(xmat, xtx, pbetamat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */
    
    SEXP betamvlsstep(const mat& xmat, const mat& xtx, const mat& wcase, const mat& hmat, 
                const colvec& y, colvec& beta, const mat& umat, double alpha,
                double taue, const colvec& zb, int nf)
    {
        BEGIN_RCPP
        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb; /* nc x 1 */
        int i; int nv = umat.n_cols;
        for(i = 0; i < nv; i++)
        {
            c += hmat.col(i) % (wcase*umat.col(i));
        }
        rmvnlschol(xmat, xtx, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

    SEXP umvstep(const mat& xmat, const mat& omega, const mat& wcase,
            const mat& wpers, const mat& hmat, const field<mat>& hmatws, 
	    const field<mat>& hwsthws, const colvec& zb,
            const colvec& beta, const colvec& y, const colvec& omegaplus,
            mat& umat, mat& mmmat, double alpha,
            double taue, const mat& L, int ns, int nc)
    {
        BEGIN_RCPP
	int nv = umat.n_cols, i;
        colvec ws(nc); colvec ytilde(nc);
        colvec es(nv), hs(nv); mat phis(nv,nv);
        colvec us_c(nv); 

	// Add in full MM term for cwithalls
	// nc x 1, c = \alpha + X\beta + zb + sum_{j=1}^{nv}[h_{j} % (W\gamma_{j})], where h_{j} is nc x 1, W is nc x ns, \gamma_{j} i ns x 1
	// this is equivalent to c_{i} = \alpha + x_{i}'\beta + zb_{i} + w_{i}'\Gamma h_{i}, for ns x 1, w_{i}, ns x nv \Gamma and nv x 1, h_{i}
	colvec cwithalls = alpha + xmat*beta + zb; /* nc x 1 */
        for(i = 0; i < nv; i++)
        {
              cwithalls 	+= hmat.col(i) % (wcase*umat.col(i));
        }
	// q, S x 1 products of S x S, omega and S x 1, omega[,q]
	mat omega_uall		= omega * umat; /* ns x nv */

        // loop over s to sample umat[s]|umat[-s],..,y
        int s, j;
        for(s = 0; s < ns; s++)
        {
            // Set entry s for data to 0
            ws = wcase.col(s); /* ws = W.case[,s] */

            // remove entry Z[,j] % wcase[,s]*u[s,j] from cwithalls
	    // remove influence of 1 x nv, u[s] from ns x nv composition of ns x ns, omega and ns x nv, umat
	    // subtract off sum_{j=1}^{q}[h_{j} % (ws\gamma_{sj})] = ws %(\sum_{j=1}^{nv}[h_{j}\gamma_{sj}] = ws % [H\gamma_{s}], where ws is column s of W and \gamma_{s} is nv x 1
	    cwithalls 	-= ( hmat * trans(umat.row(s)) ) % ws; /* [h_(1)u(s,1) + ... + h_(nv)u(s,nv)] % ws */
	    omega_uall 	-= omega.col(s) * umat.row(s); /* ns x nv */
            // sample umat[s,]
            // construct posterior mean, hs, and precision, phis
            ytilde = y - cwithalls;
            /* nv x 1 */
	    // note: omega_uall.row(s) is 1 x nv, t(omega_sbar)*omegaplus(s)
	    es 			= trans( taue*(trans(ytilde)*hmatws(s)) + omega_uall.row(s)*L );
            phis 		= taue*hwsthws(s) + omegaplus(s)*L; /* nv x nv */
            us_c.zeros(); /* nv x 1 */
            rmvnbasic(phis,es,us_c);
            umat.row(s) 	= us_c.t();
	    // puts back session s contribution from newly sampled value for 1 x nv, u(s)
            cwithalls 		+= ( hmat * trans(umat.row(s)) ) % ws;
	    omega_uall		+= omega.col(s) * umat.row(s);
        } /* end loop over s for sampling umat */
        // post-process umat to subtract out grand mean since CAR prior
        // identified only up to difference in means
        rowvec meang(nv); meang = mean(umat,0); /* dim = 0 is by column */
        mat Meang = repmat(meang,ns,1); /* ns x nv */
        umat -= Meang;

        // compute npcbt x 1, mm, resulting from mapping u to the npcbt clients
        for(j = 0; j < nv; j++)
        {
            mmmat.col(j) = wpers*umat.col(j);
        }

        END_RCPP
    } /* end function ustep to sample u[1,...,s] */


    SEXP uindmvstep(const mat& xmat, const mat& wcase,
            const mat& wpers, const mat& hmat, const field<mat>& hmatws, 
	    const field<mat>& hwsthws, const colvec& zb,
            const colvec& beta, const colvec& y, 
            mat& umat, mat& mmmat, double alpha,
            double taue, const mat& L, int ns, int nc)
    {
        BEGIN_RCPP
        int nv = umat.n_cols, i;
        colvec ws(nc); colvec ytilde(nc);
        colvec es(nv), hs(nv); mat phis(nv,nv);
        colvec us_c(nv); 

	// Add in full MM term for cwithalls
	colvec cwithalls = alpha + xmat*beta + zb; /* nc x 1 */
        for(i = 0; i < nv; i++)
        {
              cwithalls 	+= hmat.col(i) % (wcase*umat.col(i));
        }

        // loop over s to sample umat[s]|umat[-s],..,y
        int s, j;
        for(s = 0; s < ns; s++)
        {
            // Set entry s for data to 0
            ws = wcase.col(s); /* ws = W.case[,s] */

            // remove entry Z[,j] % wcase[,s]*u[s,j] from cwithalls
	    // remove influence of 1 x nv, u[s] from ns x nv composition of ns x ns, omega and ns x nv, u
	    cwithalls 	-= ( hmat * trans(umat.row(s)) ) % ws; /* [h_(1)u(1) + ... + h_(nv)u(nv)] % ws */
            // sample umat[s,]
            // construct posterior mean, hs, and precision, phis
            ytilde 		= y - cwithalls;
            /* nv x 1 */
	    es 			= taue*trans(hmatws(s))*ytilde; /* nv x 1 */
            phis 		= taue*hwsthws(s) + L; /* nv x nv */
            us_c.zeros();  /* nv x 1 */
            rmvnbasic(phis,es,us_c);
            umat.row(s) 	= us_c.t();
	    // puts back session s contribution from newly sampled value for 1 x nv, u(s)
            cwithalls 		+= ( hmat * trans(umat.row(s)) ) % ws;
        } /* end loop over s for sampling umat */

	// compute npcbt x 1, mm, resulting from mapping u to the npcbt clients
        for(j = 0; j < nv; j++)
        {
            mmmat.col(j) 	= wpers*umat.col(j);
        }

    	END_RCPP
    } /* end function ustep to sample u[1,...,s] */

    SEXP alphamvstep(const mat& xmat, const mat& wcase, const colvec& beta,
                const mat& hmat, const colvec& zb, const colvec& y,
                const mat& umat, colvec& resid, double& alpha, double taue,
                double taualph, long nc)
    {
        BEGIN_RCPP
        colvec c = xmat*beta + zb;
        int i; int nv = umat.n_cols;
        for(i = 0; i < nv; i++)
        {
            c += hmat.col(i) % (wcase*umat.col(i));
        }
        resid = y - c;
        double ea = taue*sum(resid);
        double phia = taue*double(nc) + taualph;
        double ha = ea*(1/phia);
        alpha = rnorm( 1, ha, sqrt(1/phia) )[0];
        resid -= alpha;
        END_RCPP
    } /* end function to sample intercept, alpha */
    
    SEXP alphamvlsstep(const mat& xmat, const mat& wcase, const colvec& beta,
                const mat& hmat, const colvec& zb, const colvec& y,
                const mat& umat, colvec& resid, double& alpha, double taue,
                long nc)
    {
        BEGIN_RCPP
        colvec c = xmat*beta + zb;
        int i; int nv = umat.n_cols;
        for(i = 0; i < nv; i++)
        {
            c += hmat.col(i) % (wcase*umat.col(i));
        }
        resid = y - c;
        double ea = taue*sum(resid);
        double phia = taue*double(nc);
        double ha = ea*(1/phia);
        alpha = rnorm( 1, ha, sqrt(1/phia) )[0];
        resid -= alpha;
        END_RCPP
    } /* end function to sample intercept, alpha */

    SEXP taumvstep(const mat& bmat, const colvec& resid, const mat& qmat,
            const mat& umat, const mat& V, const colvec& beta,
            mat& L, rowvec& taub, double& taue,
            double& taualph, double& taubeta, double alpha, double nu,
            double a2, double a4, double a5, double a6,
            double b2, double b4, double b5, double b6, int ns,
            int np, int nc, int nf, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* L */
        mat Vstar = inv( trans(umat)*qmat*umat + inv(V) );
        double nustar = double(ns) - double(ng) + nu;
        wishrnd(L, Vstar, nustar);
        /* taub */
        int nr = bmat.n_cols;
        a = 0.5*double(np) + a2;
        int i;
        for(i = 0; i < nr; i++)
        {
            b =  0.5*dot( bmat.col(i), bmat.col(i) ) + b2;
            taub(i) = rgamma(1, a, (1/b))[0];
        }
        /* tau1 */
        //mat bless0mat = bmat.cols(1,nr-1);
        //b = 0.5*as_scalar( sum( diagvec(bless0mat*trans(bless0mat)) ) ) + b3;
        //tau1 = rgamma(1, a, (1/b))[0];
        /* taue */
        a = 0.5*double(nc) + a4;
        b = 0.5*dot( resid, resid ) + b4;
        taue = rgamma(1, a, (1/b))[0];
        /* taualph */
        a = 0.5 + a5;
        b = 0.5*pow(alpha,2) + b5;
        taualph = rgamma(1, a, (1/b))[0];
        /* taubeta */
        a = 0.5*double(nf) + a6;
        b = 0.5*dot(beta,beta) + b6;
        taubeta = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */
    
    SEXP taumvlsstep(const mat& bmat, const colvec& resid, const mat& qmat,
            const mat& umat, const mat& V, mat& L, rowvec& taub, double& taue,
            double nu, double a2, double a4, double b2, double b4, int ns,
            int np, int nc, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* L */
        mat Vstar = inv( trans(umat)*qmat*umat + inv(V) );
        double nustar = double(ns) - double(ng) + nu;
        wishrnd(L, Vstar, nustar);
        /* taub */
        int nr = bmat.n_cols;
        a = 0.5*double(np) + a2;
        int i;
        for(i = 0; i < nr; i++)
        {
            b =  0.5*dot( bmat.col(i), bmat.col(i) ) + b2;
            taub(i) = rgamma(1, a, (1/b))[0];
        }
        /* taue */
        a = 0.5*double(nc) + a4;
        b = 0.5*dot( resid, resid ) + b4;
        taue = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */

    SEXP wishrnd(mat& L, const mat& V, double nu)
    {
        BEGIN_RCPP
        int p = L.n_rows;
        int i, j;
        mat R = chol(V);
        mat A(p,p); A.zeros();
        for(i = 0; i < p; i++)
        {
            A(i,i) = sqrt( rchisq( 1, (nu-double(i)) )[0] );
            for(j = 0; j < i; j++)
            {
                if(i  >  0)
                {
                    A(i,j) = rnorm(1,0,1)[0];
                }
            }
        }
        L = trans(R) * A * trans(A) * R;
       END_RCPP
    } /* end function rmvnqr for sampling mvn using QR decomposition */

    SEXP dic1mvcomp(const colvec& Deviance, const colvec& Alpha,
            const mat& Beta, const mat& B, const colvec& Taue,
            const mat& U, const mat& xmat, const mat& zmat, const mat& hmat,
            const mat& wcase, const colvec& y, icolvec& persons, colvec& devres,
            int np)
    {
        BEGIN_RCPP
        // compute MAP parameter estimates
        double dbar = mean(Deviance);
        double alpha = mean(Alpha);
        double taue = mean(Taue);
        rowvec beta = mean(Beta,0);
        rowvec b   = mean(B,0);
        int ns = wcase.n_cols;  int nv = U.n_cols/ns;
        rowvec u = mean(U,0); /* 1 x nv*ns */
        mat umat(ns,nv);
        int i, j;
        for(i = 0; i < nv; i++)
        {
            umat.col(i)  = trans( u.cols(i*ns,((i+1)*ns-1)) );
        }
        // compute zbhat
        int nc = xmat.n_rows;
        colvec zb(nc); zb.zeros();
        zbcomp(b, zmat, persons, zb, np);
        colvec resid = y - alpha - xmat*trans(beta) - zb;
        for(j = 0; j < nv; j++)
        {
            resid -= hmat.col(j) % (wcase*umat.col(j));
        }
        double dhat = dev(resid,taue);
        double dic = 2*dbar - dhat;
        double pd = dic - dbar;
        devres(0) = dic; devres(1) = dbar; devres(2) = dhat; devres(3) = pd;
        END_RCPP
    }

     SEXP prodHW(const mat& wcase, const mat& hmat, field<mat>& hmatws, field<mat>& hwsthws)
     {
	BEGIN_RCPP
	// schur multiply mv session effects nc x nv design matrix, hmat, by each column, ws, of wcase to produce an nc x nv matrix
	// produce quadratic product of that matrix.
	int ns	= wcase.n_cols;
	int nv = hmat.n_cols;
	int nc = hmat.n_rows; int i;
	colvec ws(ns);	mat hws(nc,nv);
	for( int s = 0; s < ns; s++ )
	{
		hws.zeros();
		ws		= wcase.col(s);
	 	for(i = 0; i < nv; i++)
               		hws.col(i) = hmat.col(i) % ws;
		hmatws(s,0)	= hws;
		hwsthws(s,0)	= hws.t() * hws;
        }
	END_RCPP
    } /* end function prodMatTwo to produce quadratic product */
