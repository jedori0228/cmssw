/** \file
 * 
 *  $Date: 2014/02/02 22:12:32 $
 *  $Revision: 1.0 $
 *
 * \author Marcello Maggi
 */


#include "DataFormats/GEMDigi/interface/GEMDigiPreReco.h"
#include <iostream>

//GEMDigiPreReco::GEMDigiPreReco (float x, float y, float ex, float ey, float corr, float tof) :
GEMDigiPreReco::GEMDigiPreReco (float x, float y, float ex, float ey, float corr, float tof, int pdgid, bool prompt) :
  x_(x),
  y_(y),
  ex_(ex),
  ey_(ey),
  corr_(corr),
  tof_(tof),
  pdgid_(pdgid),
  prompt_(prompt)
{}

GEMDigiPreReco::GEMDigiPreReco ():
  x_(0.),
  y_(0.),
  ex_(0.),
  ey_(0.),
  corr_(0.),
  tof_(-1.),
  pdgid_(0),
  prompt_(0)
{}


// Comparison
bool GEMDigiPreReco::operator == (const GEMDigiPreReco& digi) const
{
  return x_ == digi.x() and y_ == digi.y() and tof_ == digi.tof();
}


// Comparison
bool GEMDigiPreReco::operator != (const GEMDigiPreReco& digi) const
{
  return x_ != digi.x() or y_ != digi.y() or tof_ != digi.tof();
}


///Precedence operator
bool GEMDigiPreReco::operator<(const GEMDigiPreReco& digi) const
{

   if (digi.tof() == tof_){
    if(digi.x() == x_)
      return digi.y() < y_;
    else 
      return digi.x() < x_;
  } else {
    return digi.tof() < tof_;
  }
}


std::ostream & operator<<(std::ostream & o, const GEMDigiPreReco& digi)
{
//  return o << "local x=" << digi.x() << " cm y=" << digi.y()<<" cm ex=" << digi.ex() << " cm ey=" << digi.ey()<< " cm tof="<<digi.tof()<<" ns";
  return o <<  "local x=" << digi.x() << " cm y=" << digi.y()<<" cm ex=" << digi.ex() << " cm ey=" << digi.ey()<< " cm tof="<<digi.tof()<<" ns"<<" pdgID "<<digi.pdgid()<<" prompt? "<<digi.prompt();
}

void GEMDigiPreReco::print() const
{
//  std::cout << "local x=" << this->x() << " cm y=" << this->y() <<" cm tof="<<this->tof()<<" ns"<<std::endl;
  std::cout << "local x=" << this->x() << " cm y=" << this->y() <<" cm tof="<<this->tof()<<" ns"<<" pdgID "<<this->pdgid()<<" prompt? "<<this->prompt()<<std::endl;
}
