
#include "MaterialDataBase.h"

#include <stdexcept>
#include <vector>
#include <algorithm>

#include "TMaterial.h"

bool MaterialDataBase::_isInitialised = false ;

MaterialDataBase::~MaterialDataBase(){
  
  std::map<std::string,TMaterial* >::iterator it = _material_map.begin();
  std::vector<TMaterial*> deleted_objects;
  
  for( /**/; it!=_material_map.end(); ++it) 
    
    if( std::find( deleted_objects.begin(), deleted_objects.end(), (*it).second ) != deleted_objects.end() ) {
      delete (*it).second ;
      deleted_objects.push_back((*it).second) ;
    }
}

TMaterial* MaterialDataBase::getMaterial(std::string mat_name){
  
  std::map<std::string,TMaterial* >::iterator it = _material_map.find(mat_name) ;        
  
  if ( it == _material_map.end() ) { 
    MaterialDataBaseException exp;
    throw exp ; 
  } 
  else { 
    return (*it).second ; 
  }
  
}

void MaterialDataBase::initialise(){
  
  this->createMaterials(); 
  _isInitialised = true ;
  
}

void MaterialDataBase::addMaterial(TMaterial* mat, std::string name) {
  std::map<std::string, TMaterial*>::iterator it = _material_map.find(name) ; 

  MaterialDataBaseException exp;

  if ( it != _material_map.end() ) { 
    throw exp; 
  } 
  else { 
    _material_map[name] = mat  ; 
  }
}


void MaterialDataBase::createMaterials(){
  
  Double_t A, Z, density, radlen ;
  std::string name;
  
  // Beam
  A       = 14.00674 * 0.7 + 15.9994 * 0.3 ;
  Z       = 7.3 ;
  density = 1.0e-25 ; // density set to very low value
  radlen  = 1.0e25 ;  // give huge radiation length 
  name    = "beam" ;
  
  TMaterial& beam = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&beam, name);
  
  // Air
  A       = 14.00674 * 0.7 + 15.9994 * 0.3 ;
  Z       = 7.3 ;
  density = 1.205e-3 ; // g/cm^3
  radlen  = 3.42e4 ;   // cm // SJA:FIXME using the standard formular for the radiation length and the above values gives 3.06e4
  name    = "air" ;
  
  TMaterial &air = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&air, name) ;
  
  // Si
  A       = 28.09 ;
  Z       = 14.0 ;
  density = 2.33 ; // g/cm^3
  radlen  = 9.36607 ;   // cm 
  name    = "silicon";
  
  TMaterial &silicon = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&silicon, name);
  
  // C
  A       = 12.01 ;
  Z       = 6.0 ;
  density = 2.00 ; // g/cm^3
  radlen  = 21.3485 ;   // cm 
  name    = "carbon" ;
  
  TMaterial &carbon = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&carbon, name);
  
  // Aluminium
  A       = 26.9815 ;
  Z       = 13.0 ;
  density = 2.699 ; // g/cm^3
  radlen  = 24.01 ;   // cm 
  name    = "aluminium" ;
  
  TMaterial &aluminium = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.) ;
  this->addMaterial(&aluminium, name);
  
  
  // TPC Gas
  A       = 39.948*0.9+(12.011*0.2+1.00794*0.8)*0.1;
  Z       = 16.4;
  density = 0.749e-3 ;
  radlen  =  1.196e4*2;
  name    = "tpcgas" ;
  
  TMaterial &tpcgas = *new TMaterial(name.c_str(), "", A, Z, density, radlen, 0.);
  this->addMaterial(&tpcgas, name);
  
  
}

