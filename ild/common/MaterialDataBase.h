#ifndef MaterialDataBase_h
#define MaterialDataBase_h

#include <string>
#include <map>
#include <exception>

class TMaterial;

class MaterialDataBaseException: public std::exception
{
  virtual const char* what() const throw()
  {
  return "MaterialDataBaseException occurred";
  }
} ;


class MaterialDataBase 
{
  
public:
  
  static MaterialDataBase& Instance()
  {
  static MaterialDataBase singleton;
  
  if( ! _isInitialised ){
    singleton.initialise() ;
  }
  
  return singleton;
  
  }
  
  // Other non-static member functions
  
public:
  
  ~MaterialDataBase();   
  
  TMaterial* getMaterial(std::string mat_name) ;  
  
  
private:
  
  MaterialDataBase() { _material_map.clear() ;}                               // Private constructor
  
  void initialise() ;
  
  MaterialDataBase(const MaterialDataBase&) ;                 // Prevent copy-construction
  MaterialDataBase& operator=(const MaterialDataBase&) ;      // Prevent assignment
  
  void addMaterial(TMaterial* mat, std::string name); 
  void createMaterials();
  
  // private memeber variables
  std::map<std::string,TMaterial* > _material_map;
  
  static bool _isInitialised;
  
};



#endif
