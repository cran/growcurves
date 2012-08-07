/* 
 * File:   mmC.h
 * Author: savitsky
 *
 * Created on December 17, 2010, 3:31 PM
 */

#ifndef MMC_H
#define	MMC_H

#include <RcppArmadillo.h>
#include <time.h>

RcppExport SEXP DDP(SEXP yvec, SEXP Xmat, SEXP Zmat,  SEXP personsvec, SEXP Dosemat, 
           SEXP numTreat, SEXP typeTreat, SEXP OmegaList, SEXP omegapluslist,
           SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP shapealph, SEXP ratebeta, SEXP Minit);
RcppExport SEXP mmmult(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP WcaseList, SEXP MmatList, 
           SEXP OmegaList,SEXP omegaplusvecList, SEXP ngsvec, SEXP personsvec,  
           SEXP typeTreat, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP shapealph, SEXP ratebeta,
           SEXP ustrengthd);
RcppExport SEXP DPre(SEXP yvec, SEXP Xmat, SEXP Zmat,  SEXP personsvec, SEXP niterInt,
           SEXP nburnInt, SEXP nthinInt, SEXP shapealph, SEXP ratebeta);
RcppExport SEXP lgm(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP personsvec, SEXP niterInt,
           SEXP nburnInt, SEXP nthinInt);
RcppExport SEXP mmC(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson,
           SEXP Omega, SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec,
           SEXP niterInt, SEXP nburnInt, SEXP ustrengthd);
RcppExport SEXP mmCchol(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson,
           SEXP Omega, SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec,
           SEXP niterInt, SEXP nburnInt, SEXP ustrengthd);
RcppExport SEXP mmCsep(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson,
           SEXP Omega, SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec,
           SEXP niterInt, SEXP nburnInt, SEXP ustrengthd);
RcppExport SEXP mmCmiar(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Hmat, SEXP Wcase,
           SEXP Wperson, SEXP Omega, SEXP omegaplusvec, SEXP groupsvec,
           SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP nthinInt,
           SEXP ustrengthd, SEXP corsessInt, SEXP typeMM);
RcppExport SEXP mmCwb(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson,
           SEXP Omega, SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec,
           SEXP niterInt, SEXP nburnInt, SEXP ustrengthd);
RcppExport SEXP mmI(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson,
           SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP ustrengthd);
RcppExport SEXP mmIgroup(SEXP yvec, SEXP Xmat,
           SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Mmat,
           SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP ustrengthd);
RcppExport SEXP mmCplusDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Omega,
           SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec, SEXP niterInt,
           SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd, SEXP shapealph, SEXP ratebeta);
RcppExport SEXP mmIplusDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson,
           SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd,
           SEXP shapealph, SEXP ratebeta);
RcppExport SEXP mmIgroupDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Mmat,
           SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd,
           SEXP shapealph, SEXP ratebeta);
RcppExport SEXP mmCmvplusDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Hmat, SEXP Wcase, SEXP Wperson,
           SEXP Omega, SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec,
           SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd, SEXP corsessInt,
           SEXP shapealph, SEXP ratebeta, SEXP typeMM);
RcppExport SEXP growthCurve(SEXP Bmat, SEXP Alphavec, SEXP Betamat, SEXP Umat, SEXP Wperson,
           SEXP ttarmvec, SEXP personvec, SEXP TInt, SEXP maxTint, SEXP nthinInt,
           SEXP nwavesInt, SEXP modelsessions);
RcppExport SEXP predmmplusdp(SEXP ytvec, SEXP Xtmat, SEXP Ztmat, SEXP Wtcase, SEXP U,
           SEXP B, SEXP Concvec, SEXP Alphavec, SEXP Taub,
           SEXP Betamat, SEXP personstvec, SEXP rateSamp);
RcppExport SEXP preddp(SEXP ytvec, SEXP Xtmat, SEXP Ztmat,
           SEXP B, SEXP Concvec, SEXP Alphavec, SEXP Taub,
           SEXP Betamat, SEXP personstvec, SEXP rateSamp);
RcppExport SEXP predmm(SEXP ytvec, SEXP Xtmat, SEXP Ztmat, SEXP Wtcase, SEXP U,
           SEXP Alphavec, SEXP Taub,
           SEXP Betamat, SEXP personstvec, SEXP rateSamp);
RcppExport SEXP predmmmv(SEXP ytvec, SEXP Xtmat, SEXP Ztmat, SEXP Wtcase, SEXP U,
           SEXP Alphavec, SEXP Taub,
           SEXP Betamat, SEXP personstvec, SEXP rateSamp);
SEXP clusterstep(const arma::mat& xmat, const arma::mat& zmat,
            arma::field<arma::mat>& zsplit,
            arma::field<arma::colvec>& ytilsplit,
            arma::mat& pbmat, const arma::colvec& y,
            const arma::colvec& beta, double alpha, double taue,
            const arma::rowvec& taub, arma::icolvec& persons,
            arma::mat& bstarmat, arma::icolvec& s, arma::icolvec& num,
            int& M, double& conc, int np, int nr);
SEXP clustermmstep(const arma::mat& xmat, const arma::mat& zmat,
            arma::field<arma::mat>& zsplit,
            arma::field<arma::colvec>& ytilsplit, arma::mat& pbmat,
            const arma::colvec& y, const arma::mat& wcase,
            const arma::colvec& u,const arma::colvec& beta,
            double alpha, double taue, const arma::rowvec& taub,
            arma::icolvec& persons, arma::mat& bstarmat,
            arma::icolvec& s, arma::icolvec& num, int& M, double& conc,
            int np, int nr);
SEXP clustermvstep(const arma::mat& xmat, const arma::mat& zmat,
            arma::field<arma::mat>& zsplit,
            arma::field<arma::colvec>& ytilsplit, const arma::mat& hmat,
            arma::mat& pbmat,
            const arma::colvec& y, const arma::mat& wcase,
            const arma::mat& umat,const arma::colvec& beta,
            double alpha, double taue, const arma::rowvec& taub,
            arma::icolvec& persons, arma::mat& bstarmat,
            arma::icolvec& s, arma::icolvec& num, int& M, double& conc,
            int np, int nr);
SEXP concstep(double& conc, int M, int np,  double a6, double b6);
SEXP bstep(const arma::mat& xmat, const arma::mat& zmat, const arma::mat& wcase,
               const arma::colvec& y, const arma::colvec& beta,
               const arma::colvec& u, double alpha, double taue,
               const arma::rowvec& taub, arma::icolvec& persons, arma::colvec& zb,
               arma::mat& bmat, int np, int nr);
SEXP bmvstep(const arma::mat& xmat, const arma::mat& zmat, const arma::mat& hmat,
               const arma::mat& wcase, const arma::colvec& y,
               const arma::colvec& beta, const arma::mat& umat, double alpha,
               double taue, const arma::rowvec& taub, arma::icolvec& persons,
               arma::colvec& zb, arma::mat& bmat, int np, int nr);
SEXP bcholstep(const arma::mat& xmat, const arma::mat& zmat, const arma::mat& wcase,
               const arma::colvec& y, const arma::colvec& beta,
               const arma::colvec& u, double alpha, double taue,
               const arma::rowvec& taub, arma::icolvec& persons, arma::colvec& zb,
               arma::mat& bmat, int np, int nr);
SEXP blgmstep(const arma::mat& xmat, const arma::mat& zmat,
               const arma::colvec& y, const arma::colvec& beta,
               double alpha, double taue,
               const arma::rowvec& taub, arma::icolvec& persons, arma::colvec& zb,
               arma::mat& bmat, int np, int nr);
SEXP bsepstep(const arma::mat& xmat, const arma::mat& zmat, const arma::mat& wcase,
               const arma::colvec& y, const arma::colvec& beta,
               const arma::colvec& u, double alpha, double taue,
               const arma::rowvec& taub, arma::icolvec& persons, arma::colvec& zb,
               arma::mat& bmat, int np, int nr, int npcbt);
SEXP bstarstep(const arma::field<arma::mat>& zsplit,
            const arma::field<arma::colvec>& ytilsplit,
            double taue, const arma::mat& pbmat, const arma::icolvec& s,
            const arma::icolvec& num, arma::mat& bstarmat,
            arma::mat& zb, arma::mat& bmat, int np, int nr);
SEXP bpredstep(const arma::rowvec& bkvec,
            const arma::mat& ztmat, const arma::icolvec& personst,
            arma::mat& btmatk, double conck, const arma::rowvec& taubk,
            arma::colvec& zb, int npt, int np);
SEXP bmmpredstep(const arma::mat& ztmat,
            const arma::icolvec& personst, arma::mat& btmat,
            const arma::rowvec& taubk, arma::colvec& zb,
            int npt);
SEXP betastep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& y,
                 arma::colvec& beta, const arma::colvec& u, double alpha, double taue,
                 double taubeta, const arma::colvec& zb, int nc, int nf);
SEXP betalsstep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& y,
                 arma::colvec& beta, const arma::colvec& u, double alpha, double taue,
                 const arma::colvec& zb, int nc, int nf);
SEXP betamvstep(const arma::mat& xmat, const arma::mat& wcase,
                const arma::mat& hmat, const arma::colvec& y,
                arma::colvec& beta, const arma::mat& umat,
                double alpha, double taue, double taubeta,
                const arma::colvec& zb, int nf);
SEXP betamvlsstep(const arma::mat& xmat, const arma::mat& wcase,
                const arma::mat& hmat, const arma::colvec& y,
                arma::colvec& beta, const arma::mat& umat,
                double alpha, double taue, 
                const arma::colvec& zb, int nf);
SEXP betacholstep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& y,
                 arma::colvec& beta, const arma::colvec& u, double alpha, double taue,
                 double taubeta, const arma::colvec& zb, int nf);
SEXP betalscholstep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& y,
                 arma::colvec& beta, const arma::colvec& u, double alpha, double taue,
                 const arma::colvec& zb, int nf);
SEXP betalgmstep(const arma::mat& xmat, const arma::colvec& y,
                 arma::colvec& beta, double alpha, double taue,
                 const arma::colvec& zb, int nf);
SEXP betaDPstep(const arma::mat& xmat, const arma::colvec& y,
                arma::colvec& beta, double alpha, double taue,
                const arma::colvec& zb, int nf);
SEXP ustep(const arma::mat& xmat, const arma::mat& omega, const arma::mat& wcase,
               const arma::mat& wpers, const arma::colvec& beta,
               const arma::colvec& zb, const arma::colvec& y,
               const arma::colvec& omegaplus, arma::colvec& u,
               arma::colvec& mm,
               double alpha, double taue, double tauu, int ns, int nc);
SEXP umvstep(const arma::mat& xmat, const arma::mat& omega, const arma::mat& wcase,
               const arma::mat& wpers, const arma::mat& hmat,
               const arma::colvec& zb, const arma::colvec& beta,
               const arma::colvec& y, const arma::colvec& omegaplus,
               arma::mat& umat, arma::mat& mmmat, double alpha,
               double taue, const arma::mat& L, int ns, int nc);
SEXP uindmvstep(const arma::mat& xmat, const arma::mat& wcase,
               const arma::mat& wpers, const arma::mat& hmat,
               const arma::colvec& zb, const arma::colvec& beta,
               const arma::colvec& y, 
               arma::mat& umat, arma::mat& mmmat, double alpha,
               double taue, const arma::mat& L, int ns, int nc);
SEXP uwbstep(const arma::mat& xmat, const arma::mat& omega, const arma::mat& wcase,
               const arma::mat& wpers, const arma::colvec& beta,
               const arma::colvec& zb, const arma::colvec& y,
               const arma::colvec& omegaplus, arma::colvec& u,
               arma::colvec& mm, arma::icolvec& groups,
               double alpha, double taue, double tauu, int ns, int ng, int nc);
SEXP uindstep(const arma::mat& xmat, const arma::mat& wcase,
            const arma::mat& wpers, const arma::colvec& beta,
            const arma::colvec& zb, const arma::colvec& y,
            arma::colvec& u, arma::colvec& mm,  double alpha, double taue,
            double tauu, int ns);
SEXP uindetastep(const arma::mat& xmat, const arma::mat& wcase, const arma::mat& mmat,
            const arma::mat& wpers, const arma::colvec& beta,
            const arma::colvec& zb, const arma::colvec& y,
            const arma::colvec& eta,
            arma::colvec& u, arma::colvec& mm,  double alpha, double taue,
            double tauu, int ns);
SEXP etastep(const arma::mat& mmat,
            const arma::colvec& u, arma::colvec &eta,
            double tauu, double taueta);
SEXP alphastep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& beta,
                  const arma::colvec& zb, const arma::colvec& y, const arma::colvec& u,
                  arma::colvec& resid, double& alpha, double taue,
                  double taualph, long nc);
SEXP alphalsstep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& beta,
                  const arma::colvec& zb, const arma::colvec& y, const arma::colvec& u,
                  arma::colvec& resid, double& alpha, double taue,
                  long nc);
SEXP alphalgmstep(const arma::mat& xmat, const arma::colvec& beta,
                  const arma::colvec& zb, const arma::colvec& y,
                  arma::colvec& resid, double& alpha, double taue,
                  long nc);
SEXP alphamvstep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& beta,
                  const arma::mat& hmat, const arma::colvec& zb,
                  const arma::colvec& y, const arma::mat& umat,
                  arma::colvec& resid, double& alpha, double taue,
                  double taualph, long nc);
SEXP alphamvlsstep(const arma::mat& xmat, const arma::mat& wcase, const arma::colvec& beta,
                  const arma::mat& hmat, const arma::colvec& zb,
                  const arma::colvec& y, const arma::mat& umat,
                  arma::colvec& resid, double& alpha, double taue,
                  long nc);
SEXP alphaDPstep(const arma::mat& xmat, const arma::colvec& beta,
                const arma::colvec& zb, const arma::colvec& y,
                arma::colvec& resid, double& alpha, double taue,
                long nc);
SEXP taustep(const arma::mat& bmat, const arma::colvec& resid, const arma::mat& qmat,
            const arma::colvec& u, const arma::colvec& beta, double& tauu,
            arma::rowvec& taub, double& taue,
            double& taualph, double& taubeta, double alpha, double a1,
            double a2, double a4, double a5, double a6, double b1,
            double b2, double b4, double b5, double b6,
            int ns, int np, int nc, int nf, int ng);
SEXP taulsstep(const arma::mat& bmat, const arma::colvec& resid, const arma::mat& qmat,
            const arma::colvec& u, double& tauu, arma::rowvec& taub, double& taue,
            double a1, double a2, double a4, double b1, double b2, double b4, 
            int ns, int np, int nc, int ng);
SEXP taulgmstep(const arma::mat& bmat, const arma::colvec& resid, 
            arma::rowvec& taub, double& taue, double a2, double a4, 
            double b2, double b4, int np, int nc);
SEXP tausepstep(const arma::mat& bmat, const arma::colvec& resid, const arma::mat& qmat,
            const arma::colvec& u, const arma::colvec& beta, double& tauu,
            arma::rowvec& taub, double& taue,
            double& taualph, double& taubeta, double alpha, double a1,
            double a2, double a4, double a5, double a6, double b1,
            double b2, double b4, double b5, double b6,
            int ns, int npcbt, int nc, int nf, int ng);
SEXP taumvstep(const arma::mat& bmat, const arma::colvec& resid, const arma::mat& qmat,
            const arma::mat& umat, const arma::mat& V, const arma::colvec& beta,
            arma::mat& L, arma::rowvec& taub, double& taue,
            double& taualph, double& taubeta, double alpha, double nu,
            double a2, double a4, double a5, double a6,
            double b2, double b4, double b5, double b6,
            int ns, int np, int nc, int nf, int ng);
SEXP taumvlsstep(const arma::mat& bmat, const arma::colvec& resid, const arma::mat& qmat,
            const arma::mat& umat, const arma::mat& V, 
            arma::mat& L, arma::rowvec& taub, double& taue,
            double nu,double a2, double a4, double b2, double b4, 
            int ns, int np, int nc, int ng);
SEXP tauwbstep(const arma::mat& bmat, const arma::colvec& resid, const arma::mat& qmat,
            const arma::colvec& u, double& tauu, arma::rowvec& taub,
            double& taue, double a1, double a2, double a3, double a4, double b1,
            double b2, double b3, double b4, int ns, int np, int nc, int ng);
SEXP tauDPstep(const arma::mat& bstarmat, const arma::colvec& resid,
            arma::rowvec& taub, double& taue, double a2, double a3, 
            double b2, double b3, int M, int nc);
SEXP tauindstep(const arma::mat& bmat, const arma::colvec& resid,
            const arma::colvec& u, const arma::colvec& beta,
            double& tauu, arma::rowvec& taub,
            double& taue, double& taualph, double& taubeta, double alpha,
            double a1, double a2, double a4, double a5, double a6,
            double b1, double b2, double b4, double b5, double b6,
            int ns, int np, int nc, int nf);
SEXP tauetastep(const arma::mat& bmat, const arma::mat& mmat,
            const arma::colvec& resid,
            const arma::colvec& u, const arma::colvec& eta,
            const arma::colvec& beta,
            double& tauu, arma::rowvec& taub,
            double& taue, double& taualph, double& taubeta, double& taueta,
            double alpha,
            double a1, double a2, double a4, double a5, double a6, double a7,
            double b1, double b2, double b4, double b5, double b6, double b7,
            int ns, int np, int nc, int nf, int ng);
SEXP taummdpstep(const arma::mat& bstarmat, const arma::colvec& resid,
            const arma::mat& qmat, const arma::colvec& u, 
            double& tauu, arma::rowvec& taub, double& taue,
            double a1, double a2, double a4,double b1, double b2, double b4, 
            int M, int ns, int nc, int ng);
SEXP taummidpstep(const arma::mat& bstarmat, const arma::colvec& resid,
            const arma::colvec& u, 
            double& tauu, arma::rowvec& taub, double& taue,
            double a1, double a2, double a4, double b1, double b2, double b4, 
            int M, int ns, int nc);
SEXP taummigdpstep(const arma::mat& bstarmat, const arma::mat& mmat,
            const arma::colvec& resid,
            const arma::colvec& u, const arma::colvec& eta,
            double& tauu, arma::rowvec& taub, double& taue, double& taueta,
            double a1, double a2, double a4, double a7, double b1, double b2, 
            double b4, double b7, int M, int ns, int nc, int ng);
SEXP taumvdpstep(const arma::mat& bstarmat, const arma::colvec& resid,
            const arma::mat& qmat, const arma::mat& umat, const arma::mat& V, 
            arma::mat& L, arma::rowvec& taub, double& taue, double nu,
            double a2, double a4, double b2, double b4,  int M,
            int ns, int nc, int ng);
SEXP rmvnqr(const arma::mat& xmat, const arma::mat& Umat, const arma::colvec& y,
            const arma::colvec& c, arma::colvec& b, int n, int p, double taue);
SEXP rmvnlsqr(const arma::mat& xmat, const arma::colvec& y,
            const arma::colvec& c, arma::colvec& b, int n, int p, double taue);
SEXP rmvnqrclust(const arma::field<arma::mat>& zsplit, const arma::mat& Umat,
            const arma::field<arma::colvec>& ytilsplit,
            arma::colvec& b, double taue);
SEXP rmvnchol(const arma::mat& xmat, const arma::mat& Pmat, const arma::colvec& y,
            const arma::colvec& c, arma::colvec& b, int p, double taue);
SEXP rmvnlschol(const arma::mat& xmat, const arma::colvec& y,
            const arma::colvec& c, arma::colvec& b, int p, double taue);
SEXP rmvnbasic(const arma::mat& phi, const arma::colvec& e, arma::colvec& b);
SEXP rmvnmean(const arma::mat& xmat, const arma::mat& Pmat, const arma::colvec& y,
            const arma::colvec& c, const arma::colvec& m,
            arma::colvec& b, int p, double taue);
SEXP rmvncholclust(const arma::field<arma::mat>& zwedge, const arma::mat& Pmat,
            const arma::field<arma::colvec>& ytilwedge,
            arma::colvec& b, double taue);
SEXP rmvnrnd(const arma::mat& Pmat, arma::mat& B, arma::rowvec& m);
SEXP rmvncolrnd(const arma::mat& Pmat, arma::mat& B, arma::colvec& m);
double dev(const arma::colvec& resid, double taue);
SEXP dmarg(const arma::colvec& resid, double taue, arma::rowvec& devmarg);
SEXP cpo(const arma::mat& Devmarg, arma::rowvec& logcpo, double& lpml);
SEXP zbcomp(const arma::rowvec& b, const arma::mat& zmat,
            arma::icolvec& persons, arma::colvec& zb, int np);
SEXP dic1comp(const arma::colvec& Deviance, const arma::colvec& Alpha,
            const arma::mat& Beta, const arma::mat& B, 
            const arma::colvec& Taue, const arma::mat& U, const arma::mat& xmat, 
            const arma::mat& zmat, const arma::mat& wcase,const arma::colvec& y, 
            arma::icolvec& persons, arma::colvec& devres, int np);
SEXP dic1lgmcomp(const arma::colvec& Deviance, const arma::colvec& Alpha,
            const arma::mat& Beta, const arma::mat& B,
            const arma::colvec& Taue, const arma::mat& xmat,
            const arma::mat& zmat, const arma::colvec& y,
            arma::icolvec& persons, arma::colvec& devres, int np);
SEXP dic1mvcomp(const arma::colvec& Deviance, const arma::colvec& Alpha,
            const arma::mat& Beta, const arma::mat& B,
            const arma::colvec& Taue, const arma::mat& U, const arma::mat& xmat,
            const arma::mat& zmat, const arma::mat& hmat,const arma::mat& wcase,
            const arma::colvec& y,arma::icolvec& persons, arma::colvec& devres,
            int np);
SEXP dic8comp(const arma::colvec& Deviance, const arma::colvec& Alpha,
            const arma::mat& Beta, const arma::mat& B,
            const arma::colvec& Taue, const arma::mat& U, const arma::mat& xmat,
            const arma::mat& zmat, const arma::mat& wcase,const arma::colvec& y,
            arma::icolvec& persons, arma::colvec& devres, int np);
SEXP dic8dpcomp(const arma::colvec& Deviance, const arma::colvec& Alpha,
            const arma::mat& Beta, const arma::mat& B,
            const arma::colvec& Taue, const arma::mat& xmat,
            const arma::mat& zmat, const arma::colvec& y,
            arma::icolvec& persons, arma::colvec& devres, int np);
SEXP dic3comp(const arma::colvec& Deviance, const arma::mat& Devmarg, arma::colvec& devres);
double loglike(const arma::colvec& resid, double taue);
double logdens(const arma::colvec& hb, const arma::mat& phib);
unsigned long rdrawone(const arma::colvec& pr, unsigned long k);
SEXP wishrnd(arma::mat& L, const arma::mat& V, double nu);
SEXP lsqcluster(const arma::imat& S, const arma::field<arma::icolvec>& Num,
            arma::ucolvec& ordscore, arma::field<arma::ucolvec>& bigS);
SEXP ymmteststep(const arma::mat& xtmat, const arma::mat& wtmat,
        double alphak, const arma::rowvec& ukvec,
        const arma::rowvec& betakvec, const arma::colvec& zb,
        arma::colvec& ypredk);
SEXP ydpteststep(const arma::mat& xtmat, double alphak,
            const arma::rowvec& betakvec, const arma::colvec& zb,
            arma::colvec& ypredk);
SEXP ymmmvteststep(const arma::mat& xtmat, const arma::mat& wtmat,
        const arma::mat& ztmat, double alphak, const arma::rowvec& ukvec,
        const arma::rowvec& betakvec, const arma::colvec& zb,
        arma::colvec& ypredk);
SEXP buildP(const arma::field<arma::mat>& pmatsmvn, const arma::field<arma::rowvec>& taus, 
           const arma::rowvec& alphacar, const arma::rowvec& taucar, 
           const arma::field<arma::mat>& omegamats, const arma::field<arma::rowvec>& omegapluses,
           const arma::icolvec& typet, arma::mat& Ptr);
SEXP buildD(arma::field<arma::mat>& doseperson, const arma::mat& dosemat, int nr);
SEXP buildQ(const arma::mat& omegamat, const arma::rowvec& omegaplus, double alpha,
                double tau, arma::mat& qmat, int& dim);
SEXP extractE(const arma::mat& dstarmat, const arma::icolvec& typet, const arma::icolvec& numt,
                arma::field<arma::mat>& dstarmats, arma::field<arma::mat>& cstarmats, 
                arma::field<arma::colvec>& cstarvecs, int nr);
SEXP logp(const arma::mat& omegamat, const arma::rowvec& omegaplus,
            const arma::field<arma::colvec>& cstarvecs, const arma::mat& lambda, double tau, 
            double alpha, int treat, double& logf, arma::mat& qmat);
SEXP initlogp(const arma::field<arma::mat>& omegamats, const arma::field<arma::rowvec>& omegapluses,
            const arma::field<arma::colvec>& cstarvecs, const arma::mat& lambda, const arma::rowvec& taucar, 
            const arma::rowvec& alphacar, const arma::icolvec& typet, arma::colvec& logfstart, 
            arma::field<arma::mat>& qmats);
SEXP zdcomp(const arma::field<arma::mat>& doseperson, const arma::mat& dmat, 
            const arma::mat& zmat, arma::mat& thetamat, arma::icolvec& persons, 
            arma::colvec& zd);
SEXP clusterddpstep(const arma::mat& xmat, const arma::mat& zmat, 
            arma::field<arma::mat>& zsplit, arma::field<arma::colvec>& ytilsplit, 
            const arma::field<arma::mat>& doseperson, const arma::mat& Pdelt, 
            const arma::colvec& y, const arma::colvec& beta, double alpha, 
            double taue, arma::icolvec& persons, arma::mat& dstarmat, 
            arma::icolvec& s, arma::icolvec& num, int& M, 
            double& conc);
SEXP dstarstep(const arma::field<arma::mat>& zsplit, const arma::field<arma::colvec>& ytilsplit,
            double taue, const arma::mat& Pdelt, const arma::field<arma::mat>& doseperson, 
            const arma::icolvec& s, const arma::icolvec& num, arma::mat& dstarmat,
            arma::mat& zd, arma::mat& dmat, arma::mat& thetamat, int nr);
SEXP precisionddpstep(arma::field<arma::mat>& pmatsmvn, arma::field<arma::rowvec>& taus, 
           arma::mat& lambda, arma::rowvec& alphacar, arma::rowvec& taucar, arma::mat& Ptr, 
           arma::mat& Pdelt, double& taue, const arma::mat& dstarmat, 
           const arma::field<arma::mat>& omegamats, 
           const arma::field<arma::rowvec>& omegapluses, const arma::colvec& resid, 
           const arma::icolvec& typet, const arma::icolvec& numt, 
           arma::colvec& logf, arma::field<arma::mat>& qmats, int nc);
SEXP betaddpstep(const arma::mat& xmat, const arma::colvec& y,
                arma::colvec& beta, double alpha, double taue,
                const arma::colvec& zd, int nf);
SEXP alphaddpstep(const arma::mat& xmat, const arma::colvec& beta,
                const arma::colvec& zd, const arma::colvec& y,
                arma::colvec& resid, double& alpha, double taue,
                long nc);
SEXP umultsteps(arma::field<arma::colvec>& us, arma::field<arma::colvec>& ustars, 
           arma::field<arma::icolvec>& nums, arma::field<arma::ucolvec>& ss, 
           arma::rowvec& concs, arma::irowvec& Ms, arma::field<arma::colvec>& etas, 
           arma::rowvec& tauus, arma::rowvec& tauetas, const arma::field<arma::mat>& omegamats, 
           const arma::field<arma::rowvec>& omegapluses,
           const arma::field<arma::mat>& wmats, const arma::field<arma::mat>& mmats, 
           const arma::icolvec& ngs, const arma::icolvec& typet, const arma::mat& xmat, 
           const arma::colvec& y, const arma::colvec& beta, const arma::colvec& zb, 
           double alpha, double taue, arma::colvec& cmm, double ustrength);
SEXP clusterustep(const arma::colvec& y, const arma::colvec& zb, const arma::mat& xmat,
            const arma::colvec& beta, double alpha, double taue, double tauu,
            const arma::field<arma::mat>& wmats, const arma::field<arma::colvec>& us, 
            arma::colvec& ustar, arma::ucolvec& s, arma::icolvec& num, int& M, 
            double& conc, arma::colvec& ytilde, int treat, arma::colvec& cmm);
SEXP ustarstep(const arma::colvec& ytilde, const arma::mat& wtreat,
	    arma::colvec& ustar, arma::field<arma::colvec>& us, const arma::ucolvec& s, 
	    double taue, double tauu, double treat);
SEXP lsqclusteru(const arma::umat& S, arma::ucolvec& ordscore, arma::umat& bigS);
SEXP ucarstep(const arma::mat& xmat, const arma::colvec& beta, const arma::colvec& zb,
            const arma::colvec& y, const arma::field<arma::mat>& wmats, 
            arma::field<arma::colvec>& us, const arma::mat& omega, 
            const arma::rowvec& omegaplus,
            double alpha, double taue, double tauu, int treat, arma::colvec& cmm);
SEXP uistep(const arma::mat& xmat, const arma::colvec& beta, const arma::colvec& zb,
            const arma::colvec& y, const arma::field<arma::mat>& wmats, 
            arma::field<arma::colvec>& us,  
            double alpha, double taue, double tauu, int treat, arma::colvec& cmm);
SEXP uigrpstep(const arma::mat& xmat, const arma::mat& mmat, const arma::colvec& beta, 
            const arma::colvec& zb, const arma::colvec& y, const arma::colvec& eta, 
            const arma::field<arma::mat>& wmats, arma::field<arma::colvec>& us, 
            double alpha, double taue, double tauu, int treat, arma::colvec& cmm);
SEXP etagrpstep(const arma::mat& mmat, arma::colvec& eta, const arma::colvec& u,
            double tauu, double taueta);
SEXP cu(const arma::colvec& zb, const arma::mat& xmat,
            const arma::colvec& beta, double alpha, 
            const arma::field<arma::mat>& wmats, const arma::field<arma::colvec>& us,
	    arma::colvec& cmm, arma::colvec& c, int treat);
SEXP clusterbstep(const arma::mat& xmat, const arma::mat& zmat, 
            const arma::field<arma::mat>& wmats,
            const arma::field<arma::colvec>& us, arma::field<arma::mat>& zsplit,
            arma::field<arma::colvec>& ytilsplit, arma::mat& pbmat, const arma::colvec& y, 
            const arma::colvec& beta, double alpha, double taue, const arma::rowvec& taub,
            arma::icolvec& persons, arma::mat& bstarmat,
            arma::icolvec& s, arma::icolvec& num, int& M, double& conc,
            int np, int nr);
SEXP betalsmultstep(const arma::mat& xmat, const arma::field<arma::mat>& wmats, 
                const arma::field<arma::colvec>& us, const arma::colvec& y,
                arma::colvec& beta, double alpha, double taue,
                const arma::colvec& zb, int nf);
SEXP alphalsmultstep(const arma::mat& xmat, const arma::field<arma::mat>& wmats, 
                const arma::field<arma::colvec>& us, const arma::colvec& beta,
                const arma::colvec& zb, const arma::colvec& y, 
                arma::colvec& resid, double& alpha, double taue,
                long nc);
SEXP taumultstep(const arma::mat& bstarmat, const arma::colvec& resid, 
            arma::rowvec& taub, double& taue, int M, int nc);


#endif	/* MMC_H */

