#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmCchol(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Omega,
        SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec, SEXP niterInt,
        SEXP nburnInt, SEXP ustrengthd)
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
    NumericMatrix Wcr(Wcase);
    NumericMatrix Wp(Wperson);
    NumericMatrix Or(Omega);
    NumericVector opr(omegaplusvec);
    IntegerVector gr(groupsvec);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    double ustrength = as<double>(ustrengthd);
    int nkeep = niter -  nburn;

    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), ns = Or.ncol();
    int ng = gr.length(); int npcbt = Wp.nrow();

    // Compute np, number of unique clients
    // algorithm assumes cases clustered by unique client
    IntegerVector dpr = diff(pr);
    icolvec diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    /* int np = accumulate(dpr.begin,dpr.end,1); */

    // compute ng, number of session groups
    IntegerVector dgr = diff(gr);
    icolvec diffgroups(dgr.begin(),ns - 1, false);
    ng = accu(diffgroups) + 1;

    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
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
    colvec Alpha(nkeep);
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    mat Beta(nkeep,nf);
    mat B(nkeep,nr*np);
    rowvec brow(np*nr); brow.zeros();
    mat MM(nkeep,npcbt);
    mat U(nkeep,ns);
    mat Taub(nkeep,nr);
    colvec Tauu(nkeep);
    //colvec Taubeta(nkeep);
    //colvec Taualph(nkeep);
    colvec Taue(nkeep);

    // Initialize parameter values
    double tauu, taubeta, taualph, taue;
    tauu = taubeta = taualph = taue = 1;
    rowvec taub(nr); taub.ones();
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    mat    bmat = randn<mat>(np,nr)*sqrt(1/taub(0));
    colvec u = randn<colvec>(ns)*sqrt(1/tauu);
    colvec zb(nc); zb.zeros();
    mat qmat = -omega; qmat.diag() = omegaplus;/* prior precision for joint u */
    colvec resid(nc); resid.zeros(); /* regression residual */
    colvec mm(npcbt); mm.zeros(); /* client effect mapped from session eff */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    /* goodness of fit related measures */
    double deviance; colvec devres(4); colvec devres3(4);
    rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;

    // set hyperparameter values
    double a1, a2, a4, b1, b2, b4;
    a1 = b1 = ustrength; /* tauu */
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */
    //a5 = 0.01; b5 = 0.01; /* taualph */
    //a6 = 0.1; b6 = 0.1; /* taubeta */

    // Conduct posterior samples and store results
    int k, l;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) cout << "Interation: " << k << endl;
        bcholstep(xmat, zmat, wcase, y, beta, u, alpha, taue, taub,
            persons, zb, bmat, np, nr);
        // betacholstep(xmat, wcase, y, beta, u, alpha, taue, taubeta, zb, nf);
        betalscholstep(xmat, wcase, y, beta, u, alpha, taue, zb, nf);
        ustep(xmat, omega, wcase, wpers, beta, zb, y, omegaplus, u, mm, alpha, taue, tauu, ns, nc);
        // alphastep(xmat, wcase, beta, zb, y, u, resid, alpha, taue, taualph, nc);
        alphalsstep(xmat, wcase, beta, zb, y, u, resid, alpha, taue, nc);
        //taustep(bmat, resid, qmat, u, beta, tauu, taub, taue,
        //        taualph, taubeta, alpha, a1, a2, a4, a5, a6, b1, b2, 
        //        b4, b5, b6, ns, np, nc, nf, ng);
        taulsstep(bmat, resid, qmat, u, tauu, taub, taue,
                a1, a2, a4, b1, b2, b4, ns, np, nc, ng);
        if(k >= nburn)
        {
            deviance =  dev(resid, taue); /* scalar double deviance */
            dmarg(resid, taue, devmarg); /* 1 x nc vector of densities*/
            int oo = k - nburn;
            Deviance(oo) = deviance;
            Devmarg.row(oo)  = devmarg;
            Beta.row(oo) = trans(beta);
            Alpha(oo) = alpha;
            for( l = 0; l < nr; l++)
            {
                brow.cols( np*l,(np*(l+1)-1) ) = trans( bmat.col(l) );
            }
            B.row(oo) = brow;
            U.row(oo) = trans(u);
            MM.row(oo) = trans(mm);
            Tauu(oo) = tauu;
            //Taubeta(oo) = taubeta;
            //Taualph(oo) = taualph;
            Taue(oo) = taue;
            Taub.row(oo) = taub;
            Resid.row(oo) = trans(resid);
        } /* end if k > burnin, record results for return */

    } /* end MCMC loop over k */

    // DIC
    dic1comp(Deviance, Alpha, Beta, B, Taue, U, xmat, zmat, wcase,
            y, persons, devres, np); /* devres = c(dic,dbar,dhat,pd) */
    dic3comp(Deviance, Devmarg, devres3);
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  Rcpp::Named("devres") = devres,
                                  Rcpp::Named("devres3") = devres3,
                                  Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("U") = U,
                                  Rcpp::Named("MM") = MM,
                                  Rcpp::Named("Tauu") = Tauu,
                                  //Rcpp::Named("Taubeta") = Taubeta,
                                  //Rcpp::Named("Taualph") = Taualph,
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
    SEXP bcholstep(const mat& xmat, const mat& zmat, const mat& wcase,
            const colvec& y, const colvec& beta,
            const colvec& u, double alpha, double taue,
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
        // of bmat through cholesky decomposition
        mat pbmat; pbmat.eye(nr,nr);
        int i;
        for(i = 0; i < nr; i++)
        {
            pbmat.diag()(i) = taub(i);
        }
        mat zj, xj, wj, yj;
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
            xj = xmat.rows(startrow,endrow);
            wj = wcase.rows(startrow,endrow);
            yj = y.rows(startrow,endrow);

            /* sample posterior of nj block of bmat*/
            cj = alpha + xj*beta + wj*u;
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

    SEXP betacholstep(const mat& xmat, const mat& wcase, const colvec& y,
                colvec& beta, const colvec& u, double alpha, double taue,
                double taubeta, const colvec& zb, int nf)
    {
        BEGIN_RCPP
        // Compute cholesky of prior precision, pbmat, for mvn sampling
        mat pbetamat; pbetamat.eye(nf,nf); pbetamat *= taubeta;

        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb + wcase*u; /* nc x 1 */
        //rmvnqr(xmat, Ubeta, y, c, beta, nc, nf, taue);
        rmvnchol(xmat, pbetamat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */
    
    SEXP betalscholstep(const mat& xmat, const mat& wcase, const colvec& y,
                colvec& beta, const colvec& u, double alpha, double taue,
                const colvec& zb, int nf)
    {
        BEGIN_RCPP
        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb + wcase*u; /* nc x 1 */
        rmvnlschol(xmat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

    SEXP rmvnsample(const mat& phi, const colvec& h, colvec& b)
    {
        BEGIN_RCPP
        // build posterior variance and mean
        int p 		= phi.n_cols;
	colvec noise 	= randn<colvec>(p);
	b 		= solve(trimatu(chol(phi)),noise) + h;
       END_RCPP
    } /* end function rmvnbasic for drawing a single mvn sample */

    SEXP rmvnbasic(const mat& phi, const colvec& e, colvec& b)
    {
        BEGIN_RCPP
        // build posterior variance and mean
        int p 		= phi.n_cols;
        colvec h   	= solve(phi,e);
	colvec noise 	= randn<colvec>(p);
	b 		= solve(trimatu(chol(phi)),noise) + h;
       END_RCPP
    } /* end function rmvnbasic for drawing a single mvn sample */

    
     SEXP rmvnchol(const mat& xmat, const mat& Pmat, const colvec& y,
            const colvec& c, colvec& b, int p, double taue)
    {
        BEGIN_RCPP
        // S = P^-1 = (U_p' * U_p) ^-1 = U_p^-1 * (U_p^-1)' = U'U
 	// note: U' does not equal U_p^-1
	// b = U_p^-1 * z + h -> cov(b) = U_p^-1 * (U_p^-1)' = S
	colvec e 	= taue*trans(xmat)*(y-c);
	mat phi 	= symmatl(taue*trans(xmat)*xmat + Pmat);
        colvec h   	= solve(phi,e);
	colvec noise 	= randn<colvec>(p);
	b 		= solve(trimatu(chol(phi)),noise) + h; // U_p x = z -> x = solve(U_p,z)
       END_RCPP
    } /* end function rmvnchol for drawing a single mvn sample */
     
     SEXP rmvnlschol(const mat& xmat, const colvec& y,
            const colvec& c, colvec& b, int p, double taue)
    {
        BEGIN_RCPP
	colvec e 	= taue*trans(xmat)*(y-c);
	mat phi 	= symmatl( taue*trans(xmat)*xmat );
        colvec h   	= solve(phi,e);
	colvec noise 	= randn<colvec>(p);
	b 		= solve(trimatu(chol(phi)),noise) + h;
       END_RCPP
    } /* end function rmvnlschol for drawing a single mvn sample */

     SEXP rmvnmean(const mat& xmat, const mat& Pmat, const colvec& y,
            const colvec& c, const colvec& m, colvec& b, int p, double taue)
    {
        BEGIN_RCPP
        // build xbig and ytildebig
	colvec e 	= taue*trans(xmat)*(y-c) + Pmat*m;
	mat phi 	= symmatl(taue*trans(xmat)*xmat + Pmat);
        colvec h   	= solve(phi,e);
	colvec noise 	= randn<colvec>(p);
	b 		= solve(trimatu(chol(phi)),noise) + h;
       END_RCPP
    } /* end function rmvnqr for sampling mvn using QR decomposition */


    SEXP rmvnrnd(const mat& Pmat, mat& B, rowvec& m) /* by rows */
    {
        BEGIN_RCPP
        // for B(n,p), Pmat(p,p), m(p,1)
        int p 		= B.n_cols;
        int n 		= B.n_rows;
	mat U_p		= trimatu(chol(Pmat));
        int i; colvec noise(p);
        for(i = 0; i< n; i++)
        {
            noise 	= randn<colvec>(p);
	    B.row(i)	= trans(solve(U_p,noise)) + m;
        }
       END_RCPP
    } /* end function rmvnrd for cholesky decomp sampling of mvn */

    SEXP rmvncolrnd(const mat& Pmat, mat& B, colvec& m) /* by cols */
    {
        BEGIN_RCPP
        // for B(n,p), Pmat(n,n), m(n,1)
        //int p 	= B.n_cols;
        int n 		= B.n_rows;
	mat U_p		= trimatu(chol(Pmat));
        int k; colvec noise(n);
        for(k = 0; k < n; k++)
        {
            noise 		= randn<colvec>(n);
            B.col(k) 		= solve(U_p,noise)  + m;
        }
       END_RCPP
    } /* end function rmvnrd for cholesky decomp sampling of mvn */
