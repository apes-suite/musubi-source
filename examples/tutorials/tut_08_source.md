title: Source Terms
@warning WORK IN PROGRESS @endwarning

Navigate: [&larr; Abort Criteria](tut_07_convergence.md)
| [Overview](index.md)
| [Multi-Level Simulations &rarr;](tut_09_mus_multilevel.md)

# Source Terms

Force as a source term.

An example testcase can be found in `examples/tutorials/tutorial_cases/tutorial_PIP_Force`.

To define a force that affects the whole domain one uses the `glob_source` table.

```lua
--! [Source]
glob_source = {
  force = {press_grad, 0.0, 0.0},
  force_order = 2
}
--! [Source]
```

The force `press_grad` has to be defined inside the `variable` table.

```lua
variable = {
  {
    name = 'press_grad',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = press_grad
  },
```

@note: `st_fun = press_grad` refers to a lua-function defined somewhere else in
the lua file.

Next chapter: [Multi-Level Simulations &rarr;](tut_09_mus_multilevel.md)
