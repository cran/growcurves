#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmC(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Omega,
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
    int nc = Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), ns = Or.ncol();
    int ng = gr.length(); int npcbt = Wp.nrow();

    // Compute np, number of unique clients
    // algorithm assumes cases clustered by unique client
    IntegerVector dpr = diff(pr);
    Col<int32_t> diffpersons(dpr.begin(), nc - 1, false);
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
    Col<int32_t> groups(gr.begin(), ns, false);
    Col<int32_t> persons(pr.begin(), nc, false);

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
    mat B(nkeep,np*nr);
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
    //colvec d(ns); d.zeros(); /* prior mean for u */
    mat qmat = -omega; qmat.diag() = omegaplus;/* prior precision for joint u */
    colvec resid(nc); resid.zeros(); /* regression residual */
    colvec mm(npcbt); mm.zeros(); /* client effect mapped from session eff */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    double deviance; colvec devres(4);
    rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;

    // pre-compute quadratic product of design matrix for use in sampling
    colvec wstws(ns); dotMMvecs(wcase, wstws); /* sequential, by effect, sampling of MM effects */

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
        bstep(xmat, zmat, wcase, y, beta, u, alpha, taue, taub,
            persons, zb, bmat, np, nr);
        //betastep(xmat, wcase, y, beta, u, alpha, taue, taubeta, zb, nc, nf);
        betalsstep(xmat, wcase, y, beta, u, alpha, taue, zb, nc, nf);
        ustep(xmat, omega, wcase, wstws, wpers, beta, zb, y, omegaplus, u, mm, alpha,
                taue, tauu, ns, nc);
        //alphastep(xmat, wcase, beta, zb, y, u, resid, alpha, taue, taualph, nc);
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
    SEXP bstep(const mat& xmat, const mat& zmat, const mat& wcase,
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
        // of bmat through QR decomposition
        mat pbmat; pbmat.eye(nr,nr);
        int i;
        for(i = 0; i < nr; i++)
        {
            pbmat.diag()(i) = taub(i);
        }
        mat Ub = chol(pbmat); /* returns upper triangler, U: t(U)*U = pbmat */
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
            rmvnqr(zj, Ub, yj, cj, bj, nj, nr, taue);
            bmat.row(persrow) = trans(bj);

            // compute zbj (jth piece) of zb = [z1*b1vec,...,znp*bnpvec]
            zb.rows(startrow,endrow) = zj*bj;

            //re-set start positions
            startrow += nj;
            persrow++; /*bmat has np rows, increment once for each iteration*/
        } /* end loop j for sampling bj*/
        END_RCPP
    } /* end function bstep for sampling bmat and zb */

    SEXP betastep(const mat& xmat, const mat& wcase, const colvec& y,
                colvec& beta, const colvec& u, double alpha, double taue,
                double taubeta, const colvec& zb, int nc, int nf)
    {
        BEGIN_RCPP
        // Compute cholesky of prior precision, pbmat, for mvn sampling
        // of bmat through QR decomposition
        mat pbetamat; pbetamat.eye(nf,nf); pbetamat *= taubeta;
        mat Ubeta = chol(pbetamat); /* R: t(U)*U = pbetamat */

        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb + wcase*u; /* nc x 1 */
        rmvnqr(xmat, Ubeta, y, c, beta, nc, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */
    
    SEXP betalsstep(const mat& xmat, const mat& wcase, const colvec& y,
                colvec& beta, const colvec& u, double alpha, double taue,
                const colvec& zb, int nc, int nf)
    {
        BEGIN_RCPP
        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb + wcase*u; /* nc x 1 */
        rmvnlsqr(xmat, y, c, beta, nc, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

    SEXP ustep(const mat& xmat, const mat& omega, const mat& wcase,
            const colvec& wstws, const mat& wpers, const colvec& beta, const colvec& zb,
            const colvec& y, const colvec& omegaplus, colvec& u, 
            colvec& mm,  double alpha, double taue,
            double tauu, int ns, int nc)
    {
        BEGIN_RCPP
        // build portion of offset constant, c, not dependent on u
        colvec cwithalls = alpha + xmat*beta + zb + wcase*u; 
	colvec omega_uall = omega*u;
        // set up structures that will set column s of wcase = 0 and
        // entry s for row s of omega = 0 for sampling u[s]
        colvec ws(nc); colvec ytilde;
        double es, hs, phis;

        // loop over s to sample u[s]|u[-s],..,y
        int s;
        for(s = 0; s < ns; s++)
        {
            ws = wcase.col(s); /* w.s = W.case[,s] */
            // remove entry W.case[,-s]*u[-s] from cwithalls
            cwithalls 	-= ws*u(s); /* removes session s MM contribution */
	    omega_uall	-= omega.col(s)*u(s); 	/* S x 1 vector of omega * u with effect of session s, u(s), removed */
            // sample u[s]
            // construct posterior mean, hs, and precision, phis
            ytilde 	= y - cwithalls;
            es 		= taue*dot(ytilde,ws) + tauu*omega_uall(s); /* omega_uall(s) is the product of row s of omega with u, after removing effect, u(s) */
            phis 	= taue*wstws(s) + tauu*omegaplus(s);
            hs 		= es*(1/phis);
            u(s) 	= rnorm( 1, hs, sqrt(1/phis) )[0];
	    // puts back session s contribution from newly sampled value for u(s)
            cwithalls 	+= ws*u(s);
	    omega_uall	+= omega.col(s)*u(s);
        } /* end loop over s for sampling u */
        // post-process u to subtract out group means since CAR prior
        // identified only up to difference in means
        u -= mean(u);

        // determine number of members per group for the ns sessions.
        // assumes groups arranged in continguous, ordered, disjoint blocks
        //icolvec positions(ng); positions.zeros();
       // icolvec::iterator aa = groups.begin();
       // icolvec::iterator bb = groups.end();
       // for(icolvec::iterator i = aa; i != bb ; i++)
       // {
       //    positions(*i-1) += 1;
       // }
        // compute mean of u for each group
        //int g = 0; int startrow = 0; double meang;
        //for(g = 0; g < ng; g++)
        //{
           // int nsg = positions(g);
          //  int endrow = startrow + nsg - 1;
           // meang = mean(u.rows(startrow,endrow));
          //  u.rows(startrow,endrow) -= meang;
          //  startrow += nsg;
        //} /* end loop g over nsg groups to subtract out group mean from u */


        // compute npcbt x 1, mm, resulting from mapping u to the npcbt clients
        mm = wpers*u;

        END_RCPP
    } /* end function ustep to sample u[1,...,s] */

    SEXP alphastep(const mat& xmat, const mat& wcase, const colvec& beta,
                const colvec& zb, const colvec& y, const colvec& u,
                colvec& resid, double& alpha, double taue,
                double taualph, long nc)
    {
        BEGIN_RCPP
        colvec c = xmat*beta + wcase*u + zb;
        resid = y - c;
        double ea = taue*sum(resid);
        double phia = taue*double(nc) + taualph;
        double ha = ea*(1/phia);
        alpha = rnorm( 1, ha, sqrt(1/phia) )[0];
        resid -= alpha;
        END_RCPP
    } /* end function to sample intercept, alpha */
    
    SEXP alphalsstep(const mat& xmat, const mat& wcase, const colvec& beta,
                const colvec& zb, const colvec& y, const colvec& u,
                colvec& resid, double& alpha, double taue,
                long nc)
    {
        BEGIN_RCPP
        colvec c = xmat*beta + wcase*u + zb;
        resid = y - c;
        double ea = taue*sum(resid);
        double phia = taue*double(nc);
        double ha = ea*(1/phia);
        alpha = rnorm( 1, ha, sqrt(1/phia) )[0];
        resid -= alpha;
        END_RCPP
    } /* end function to sample intercept, alpha */

    SEXP taustep(const mat& bmat, const colvec& resid, const mat& qmat,
            const colvec& u, const colvec& beta,
            double& tauu, rowvec& taub, double& taue,
            double& taualph, double& taubeta, double alpha, double a1,
            double a2, double a4, double a5, double a6, double b1,
            double b2, double b4, double b5, double b6, int ns,
            int np, int nc, int nf, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* tauu */
        a = 0.5*(double(ns) - double(ng)) + a1;
        b = 0.5*( as_scalar(trans(u)*qmat*u) ) + b1;
        // b = 0.5* dot( omegaplus, square(u-d) ) + b1;
        tauu = rgamma(1, a, (1/b))[0];
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
       // a = 0.5*double(np) + a3;
       // mat bless0mat = bmat.cols(1,nr-1);
       // b = 0.5*as_scalar( sum( diagvec(bless0mat*trans(bless0mat)) ) ) + b3;
        // tau1 = rgamma(1, a, (1/b))[0];
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

    SEXP taulsstep(const mat& bmat, const colvec& resid, const mat& qmat,
            const colvec& u,  double& tauu, rowvec& taub, double& taue,
            double a1, double a2, double a4, double b1, double b2, double b4, 
            int ns, int np, int nc, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* tauu */
        a = 0.5*(double(ns) - double(ng)) + a1;
        b = 0.5*( as_scalar(trans(u)*qmat*u) ) + b1;
        tauu = rgamma(1, a, (1/b))[0];
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
    } /* end function taulsstep to sample precision parameters */
    
    double dev(const colvec& resid, double taue)
    {
        // NumericVector r(resid);
        // devvec is an nc x 1 vector
        NumericVector r = wrap(resid);
        NumericVector devvec  = dnorm( r, 0.0, sqrt(1/taue), true ); /* logD */
        double deviance = accumulate(devvec.begin(),devvec.end(),0.0);
        deviance *= -2;
        return(deviance);
    }

    SEXP zbcomp(const rowvec& b, const mat& zmat,
            icolvec& persons, colvec& zb, int np)
    {
        BEGIN_RCPP
        // compute number of repeated measures per person
        // positions is np x 1 with each cell the number of measures
        icolvec positions(np); positions.zeros();
        icolvec::iterator aa = persons.begin();
        icolvec::iterator bb = persons.end();
        for(icolvec::iterator i = aa; i != bb ; i++)
        {
           positions(*i-1) += 1;
        }
        // compute zb, nc x 1
        int startrow = 0;
        int nr = zmat.n_cols;
        int i, j, nj, endrow; mat zj;
        colvec bj(nr);
        for(j=0; j < np; j++)
        {
            /* extract by-person matrix views from data*/
            nj = positions(j);
            endrow = startrow + nj - 1;
            zj = zmat.rows(startrow,endrow);
            // compute zbj (jth piece) of zb = [z1*b1vec,...,znp*bnpvec]
            for(i = 0; i < nr; i++) /* pick jth component of each 0 < i < nr */
            {
                bj(i) = b(i*np + j);
            }
            zb.rows(startrow,endrow) = zj*bj;
            //re-set start positions
            startrow += nj;
        } /* end loop j for sampling bj*/
        END_RCPP
    } /* end function to compute zb*/


    SEXP dic1comp(const colvec& Deviance, const colvec& Alpha,
            const mat& Beta, const mat& B, const colvec& Taue,
            const mat& U, const mat& xmat, const mat& zmat, const mat& wcase,
            const colvec& y, icolvec& persons, colvec& devres, int np)
    {
        BEGIN_RCPP
        // compute MAP parameter estimates
        double dbar = mean(Deviance);
        double alpha = mean(Alpha);
        double taue = mean(Taue);
        rowvec beta = mean(Beta,0);
        rowvec b = mean(B,0);
        rowvec u = mean(U,0);
        // compute zbhat
        int nc = xmat.n_rows;
        colvec zb(nc); zb.zeros();
        zbcomp(b, zmat, persons, zb, np);
        colvec resid = y - alpha - xmat*trans(beta) - wcase*trans(u) - zb;
        double dhat = dev(resid,taue);
        double dic = 2*dbar - dhat;
        double pd = dic - dbar;
        devres(0) = dic; devres(1) = dbar; devres(2) = dhat; devres(3) = pd;
        END_RCPP
    }


    SEXP rmvnqr(const mat& xmat, const mat& Umat, const colvec& y,
            const colvec& c, colvec& b, int n, int p, double taue)
    {
        BEGIN_RCPP
        // build xbig and ytildebig
        double errinv = sqrt(taue);
        mat xbig(n+p,p); xbig.zeros();
        colvec ytildebig(n+p); ytildebig.zeros();
        colvec noisevec = randn<colvec>(p);
        xbig.rows(0,n-1) = errinv*xmat;
        xbig.rows(n,n+p-1) = Umat;
        ytildebig.rows(0,n-1) = errinv*(y - c);
        // conduct QR backsolving for bcvec
        mat Q, R;
        qr(Q,R,xbig);
        // Translate from Matlab method where
        // for mxn, X where m > n, Q is m x m and R is m x n
        // The m-n cols of Q are basis for OC of X and m-n rows of R are O
        Q = Q.cols(0,p-1);
	R = R.rows(0,p-1);
        colvec q = trans(Q)*ytildebig;
        q += noisevec;
        // return b
        b = solve(R,q);
        END_RCPP
    } /* end function rmvnqr for sampling mvn using QR decomposition */
    
    SEXP rmvnlsqr(const mat& xmat, const colvec& y,
            const colvec& c, colvec& b, int n, int p, double taue)
    {
        BEGIN_RCPP
        // build xbig and ytildebig
        double errinv = sqrt(taue);
        mat xls(n,p); xls.zeros();
        colvec ytildels(n); ytildels.zeros();
        colvec noisevec = randn<colvec>(p);
        xls.rows(0,n-1) = errinv*xmat;
        ytildels.rows(0,n-1) = errinv*(y - c);
        // conduct QR backsolving for bcvec
        mat Q, R;
        qr(Q,R,xls);
        // Translate from Matlab method where
        // for mxn, X where m > n, Q is m x m and R is m x n
        // The m-n cols of Q are basis for OC of X and m-n rows of R are O
        Q = Q.cols(0,p-1);
	R = R.rows(0,p-1);
        colvec q = trans(Q)*ytildels;
        q += noisevec;
        // return b
        b = solve(R,q);
        END_RCPP
    } /* end rmvnlsqr for sampling mvn under least squares using QR decomposition */

    // double s2 = std::inner_product(y.begin(),y.end(),y.begin(),double());

    SEXP dotMMvecs(const mat& wcase, colvec& wstws)
    {
	BEGIN_RCPP
	// split MM weight matrix by MM effect and compute inner product
	int ns	= wcase.n_cols;
	colvec ws(ns);	
	for(int s = 0; s < ns; s++)
	{
		ws		= wcase.col(s);
		wstws(s)	= dot(ws,ws);
	}
	END_RCPP
    } /* end function prodMatTwo to produce quadratic product */
