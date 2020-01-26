import bpy

# Create new material
material_name = 'Test'
material = bpy.data.materials.new(name=material_name)
material.use_nodes = True

nodes = material.node_tree.nodes
nodes.clear()

absorption = nodes.new(type='ShaderNodeVolumeAbsorption')
scattering = nodes.new(type='ShaderNodeVolumeScattering')
output_material = nodes.new(type='ShaderNodeOutputMaterial')
output_world = nodes.new(type='ShaderNodeOutputMaterial')

# Custom nodes

shader RefractiveIndex (float n = 1.0,  output float Value = 1.0) {
    Value = n;
}
