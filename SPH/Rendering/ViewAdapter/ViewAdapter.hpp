//
//  ViewAdapter.hpp
//  Lattice Boltzman
//
//  Created by Charlie Close on 08/12/2024.
//

#ifndef VIEW_ADAPTER_HPP
#define VIEW_ADAPTER_HPP

#include <MetalKit/MetalKit.hpp>
#import "Camera.hpp"


namespace Explorer {
class ViewAdapter {
public:
  ViewAdapter(Camera* camera);
  virtual MTK::View* get(CGRect frame);
  virtual void printDebug();
  static Camera* _camera;
};
}; // namespace Explorer

#endif // VIEW_ADAPTER_HPP
