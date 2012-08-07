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

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    srand ( time(NULL) ); /* arma */

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
        bmvstep(xmat, zmat, hmat, wcase, y, beta, umat, alpha, taue, taub,
            persons, zb, bmat, np, nr);
        betamvlsstep(xmat, wcase, hmat, y, beta, umat, alpha, taue, 
                    zb, nf);
	if( typemm == 0)
	{
		uindmvstep(xmat, wcase, wpers, hmat, zb, beta, y, umat,
                	mmmat, alpha, taue, L, ns, nc);
	}else{ /* typemm == 1, MM-MCAR */
		umvstep(xmat, omega, wcase, wpers, hmat, zb, beta, y, omegaplus, umat,
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
    SEXP bmvstep(const mat& xmat, const mat& zmat, const mat& hmat,
            const mat& wcase,const colvec& y, const colvec& beta,
            const mat& umat, double alpha, double taue,
            const rowvec& taub, icolvec& persons, colvec& zb,
            mat& bmat, int np, int nr)
    {
        BEGIN_RCPP
        // compute number of cases for each person and store results in
        // 'positions' (np x 1)
        icolvec positions(np); positions.zeros();
        icolvec::iterator aa = persons.begin();
        icolvec::iterator bb = persons.end();
        for(icolvec::iterator i = aa; i != bb ; i++)
        {
           positions(*i-1) += 1;
        }

        // Compute cholesky of prior precision, pbmat, for mvn sampling
        // of bmat through QR decomposition
        mat pbmat; pbmat.eye(nr,nr);
        int i;
        for(i = 0; i < nr; i++)
        {
            pbmat.diag()(i) = taub(i);
        }
        mat hj, zj, xj, wj, yj;
        colvec cj, bj(nr);

        // sample bmat, independently, by person
        int startrow = 0; int persrow = 0;
        int j, nj, endrow;
        for(j=0; j < np; j++)
        {
            /* extract by-person matrix views from data*/
            nj = positions(j);
            endrow = startrow + nj - 1;
            zj = zmat.rows(startrow,endrow);
            hj = hmat.rows(startrow,endrow);
            xj = xmat.rows(startrow,endrow);
            wj = wcase.rows(startrow,endrow);
            yj = y.rows(startrow,endrow);

            /* sample posterior of nj block of bmat*/
            cj = alpha + xj*beta;
            int nv = umat.n_cols;
            for(i = 0; i < nv; i++)
            {
                cj += hj.col(i) % (wj*umat.col(i));
            }
            bj.zeros();
            //rmvnqr(zj, Ub, yj, cj, bj, nj, nr, taue);
            rmvnchol(zj, pbmat, yj, cj, bj, nr, taue);
            bmat.row(persrow) = trans(bj);

            // compute zbj (jth piece) of zb = [z1*b1vec,...,znp*bnpvec]
            zb.rows(startrow,endrow) = zj*bj;

            //re-set start positions
            startrow += nj;
            persrow++; /*bmat has np rows, increment once for each iteration*/
        } /* end loop j for sampling bj*/
        END_RCPP
    } /* end function bstep for sampling bmat and zb */

    SEXP betamvstep(const mat& xmat, const mat& wcase, const mat& hmat,
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
        //rmvnqr(xmat, Ubeta, y, c, beta, nc, nf, taue);
        rmvnchol(xmat, pbetamat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */
    
    SEXP betamvlsstep(const mat& xmat, const mat& wcase, const mat& hmat,
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
        rmvnlschol(xmat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

    SEXP umvstep(const mat& xmat, const mat& omega, const mat& wcase,
            const mat& wpers, const mat& hmat, const colvec& zb,
            const colvec& beta, const colvec& y, const colvec& omegaplus,
            mat& umat, mat& mmmat, double alpha,
            double taue, const mat& L, int ns, int nc)
    {
        BEGIN_RCPP
        // build portion of offset constant, c, not dependent on u
        colvec cminuss = alpha + xmat*beta + zb;
        // set up structures that will set column s of wcase = 0 and
        // entry s for row s of omega = 0 for sampling u[s]
        mat wcases; rowvec omegarow(ns); colvec ws(nc); colvec c(nc);
        int nv = umat.n_cols;
        colvec es(nv), hs(nv); mat phis(nv,nv);
        rowvec umatsbar(nv); mat us(1,nv); rowvec ths(nv); mat hws(nc,nv);

        // loop over s to sample umat[s]|umat[-s],..,y
        int s, i, j;
        for(s = 0; s < ns; s++)
        {
            // Set entry s for data to 0
            umat.row(s).zeros();
            omegarow = omega.row(s);
            omegarow(s) = 0;
            ws = wcase.col(s); /* w.s = W.case[,s] */
            wcases = wcase;
            wcases.col(s).zeros(); /* set W.case[,s] = 0*/
            // put entry Z[,j]*W.case[,-s]%*%U[-s,j] into c
            /* adds in nothing for session s */
            c.zeros(); /* re-zero contribution of CAR for each iteration */
            for(i = 0; i < nv; i++)
            {
                c += hmat.col(i) % (wcases*umat.col(i));
                hws.col(i) = hmat.col(i) % ws;
                /* mean umat(0), ..., umat(nv) */
                umatsbar(i) = dot( omegarow, umat.col(i) )/omegaplus(s);
            }
            c += cminuss;
            // sample umat[s,]
            // construct posterior mean, hs, and precision, phis
            colvec ytilde = y - c;
            /* nv x 1 */
	    es = trans( taue*(trans(ytilde)*hws) + omegaplus(s)*umatsbar*L );
            phis = taue*trans(hws)*hws + omegaplus(s)*L; /* nv x nv */
            hs = inv(phis)*es; /* nv x 1*/
            us = umat.row(s);
            ths = trans(hs);
            rmvnrnd(phis, us, ths);
            umat.row(s) = us.row(0);
        } /* end loop over s for sampling umat */
        // post-process umat to subtract out grand mean since CAR prior
        // identified only up to difference in means
        rowvec meang(nv); meang = mean(umat,0); /* dim = 0 is by column */
        mat Meang = repmat(meang,ns,1); /* ns x nv */
        umat -= Meang;

        // determine number of members per group for the ns sessions.
        // assumes groups arranged in continguous, ordered, disjoint blocks
        //icolvec positions(ng); positions.zeros();
        //icolvec::iterator aa = groups.begin();
        //icolvec::iterator bb = groups.end();
        //for(icolvec::iterator i = aa; i != bb ; i++)
        //{
        //   positions(*i-1) += 1;
       // }
        // compute mean of u for each group
        //int g; int startrow = 0; rowvec meang; mat Meang;
       // for(g = 0; g < ng; g++)
       // {
        //    int nsg = positions(g);
         //   int endrow = startrow + nsg - 1;
        //   meang = mean(umat.rows(startrow,endrow),0);
         //   Meang = repmat(meang,nsg,1);
         //   umat.rows(startrow,endrow) -= Meang;
          //  startrow += nsg;
       // } /* end loop g over nsg groups to subtract out group mean from u */

        // re-set d with final values for u in later posterior draw for tauu
        // mat omegas = omega; omegas.diag().zeros();
        // d = (omegas*u) / omegaplus;

        // compute npcbt x 1, mm, resulting from mapping u to the npcbt clients
        for(j = 0; j < nv; j++)
        {
            mmmat.col(j) = wpers*umat.col(j);
        }

        END_RCPP
    } /* end function ustep to sample u[1,...,s] */

    SEXP uindmvstep(const mat& xmat, const mat& wcase,
            const mat& wpers, const mat& hmat, const colvec& zb,
            const colvec& beta, const colvec& y, 
            mat& umat, mat& mmmat, double alpha,
            double taue, const mat& L, int ns, int nc)
    {
        BEGIN_RCPP
        // build portion of offset constant, c, not dependent on u
        colvec cminuss = alpha + xmat*beta + zb;
        // set up structures that will set column s of wcase = 0 and
        // entry s for row s of omega = 0 for sampling u[s]
        mat wcases; colvec ws(nc), c(nc);
        int nv = umat.n_cols;
        colvec es(nv), hs(nv); mat phis(nv,nv);
        mat us(1,nv); rowvec ths(nv); mat hws(nc,nv);

        // loop over s to sample umat[s]|umat[-s],..,y
        int s, i, j;
        for(s = 0; s < ns; s++)
        {
            // Set entry s for data to 0
            umat.row(s).zeros();
            ws = wcase.col(s); /* w.s = W.case[,s] */
            wcases = wcase;
            wcases.col(s).zeros(); /* set W.case[,s] = 0*/
            // put entry Z[,j]*W.case[,-s]%*%U[-s,j] into c
            /* adds in nothing for session s */
            c.zeros(); /* re-zero contribution of terms -s for each iteration */
            for(i = 0; i < nv; i++)
            {
                c += hmat.col(i) % (wcases*umat.col(i));
                hws.col(i) = hmat.col(i) % ws;
            }
            c += cminuss;
            // sample umat[s,]
            // construct posterior mean, hs, and precision, phis
            colvec ytilde = y - c;
            /* nv x 1 */
	    es = taue*trans(hws)*ytilde; /* nv x 1 */
            phis = taue*trans(hws)*hws + L; /* nv x nv */
            hs = inv(phis)*es; /* nv x 1*/
            us = umat.row(s);
            ths = trans(hs);
            rmvnrnd(phis, us, ths);
            umat.row(s) = us.row(0);
        } /* end loop over s for sampling umat */

	// compute npcbt x 1, mm, resulting from mapping u to the npcbt clients
        for(j = 0; j < nv; j++)
        {
            mmmat.col(j) = wpers*umat.col(j);
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
