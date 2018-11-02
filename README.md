# CS

Coordinate systems library for planetary N-body simulation.

## Usage

### Example

Convert Orbital Elements to Cartesian Coordinates

```cpp
// Arrange containers to store results
vector<valarray<double>> pos(3, valarray<double>(3));
vector<valarray<double>> vel(3, valarray<double>(3));

// Convert orbital elements to cartesian coordinates
cs::ConvertOelemToCartesian(mass, oelems, &pos, &vel);
```

Convert Cartesian Coordinates to Orbital Elements

``` cpp
auto oelems = cs::ConvertCartesianToOelem(mass, pos, vel);
```

