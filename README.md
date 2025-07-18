## Overview

This library implements concepts from algebraic topology:

- **Simplicial complexes**: Collections of simplices forming geometric structures  
- **Chain complexes**: Sequences of abelian groups connected by boundary operators  
- **Free modules**: Algebraic structures generalizing vector spaces  

## Design

### 1. Proper Encapsulation

- Private member variables with trailing underscore convention (`value_`, `generators_`, etc.)
- Public getter/setter methods for controlled access
- Const-correctness throughout the codebase

### 2. Generic Programming Excellence

- Template specializations for edge cases (`ZMod<0>`, `Kompleks<S,0,p>`)
- Iterator support with proper STL-compatible interfaces

### 3. Modern C++ Best Practices

- Rule of Five compliance (constructors, destructor, copy/move operations)
- Exception safety with proper error handling

### 4. Enhanced Functionality

- Comprehensive validation of operations and indices
- STL-compatible iterators for range-based operations
- Utility methods for common mathematical operations


(the boundary operations are still work in progress)
