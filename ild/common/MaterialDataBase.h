#ifndef MaterialDataBase_h
#define MaterialDataBase_h

/** MaterialDataBase: Class to hold and manage collection of materials  
 *
 * @author S.Aplin DESY
 */

#include <string>
#include <map>
#include <exception>

class TMaterial;

class MaterialDataBaseException: public std::exception {
  virtual const char* what() const throw() {
    return "MaterialDataBaseException occurred";
  }
} ;


class MaterialDataBase {
  
public:
  
  /** Accessor Method */
  static MaterialDataBase& Instance() {
    
    static MaterialDataBase singleton;
    
    if( ! _isInitialised ){
      singleton.initialise() ;
    }
    
    return singleton;
    
  }
  
  // Other non-static member functions
  
public:
  
  /** Destructor */
  ~MaterialDataBase();   
  
  /** Get Material via name */
  TMaterial* getMaterial(std::string mat_name) ;  
  
  
private:
  
  MaterialDataBase() { _material_map.clear() ;}                               // Private constructor
  
  void initialise() ;
  
  MaterialDataBase(const MaterialDataBase&) ;                 // Prevent copy-construction
  MaterialDataBase& operator=(const MaterialDataBase&) ;      // Prevent assignment
  
  void addMaterial(TMaterial* mat, std::string name); 
  void createMaterials();
  
  // private member variables
  std::map<std::string,TMaterial* > _material_map;
  
  static bool _isInitialised;
  
};



#endif
