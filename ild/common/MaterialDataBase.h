#ifndef MaterialDataBase_h
#define MaterialDataBase_h

#include <string>
#include <map>
#include <exception>

class TMaterial;

class MaterialDataBaseexception: public std::exception
{
  virtual const char* what() const throw()
  {
    return "MaterialDataBaseexception occurred";
    }
} ;

class MaterialDataBase 
{
  public:
  static MaterialDataBase& Instance()
  {
    static MaterialDataBase singleton;
    return singleton;
  }
  
    // Other non-static member functions
  
 public:
  
  ~MaterialDataBase();   

  TMaterial* getMaterial(std::string mat_name) ;  
  void initialise() ;
  
 private:
 MaterialDataBase(): _isInitialised(false) { _material_map.clear() ;} ;                               // Private constructor
  MaterialDataBase(const MaterialDataBase&) ;                 // Prevent copy-construction
  MaterialDataBase& operator=(const MaterialDataBase&) ;      // Prevent assignment
  
  void addMaterial(TMaterial* mat, std::string name); 
  void createMaterials();
  
  // private memeber variables
  std::map<std::string,TMaterial* > _material_map;
  
  bool _isInitialised;
  
};

#endif
