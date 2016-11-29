#include "fwdmodel_sine_lib.h"
#include "fwdmodel_sine.h"

extern "C" {
  int get_num_models()
  {
    return 1;
  }

  const char *get_model_name(int index)
  {
    switch(index) {
      case 0:
        return "sine";
        break;
      default:
        return NULL;
    }
  }

  NewInstanceFptr get_new_instance_func(const char *name)
  {
    if (string(name) == "sine") {
      return SineFwdModel::NewInstance;
    }
    else {
      return NULL;
    }
  }
}

