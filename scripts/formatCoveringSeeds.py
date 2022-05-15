
def formatCoveringSeeds():
    #hide/reveal
    cmd.hide("all")
    cmd.show("sticks","peptide"+" and bb.")
    cmd.show("sticks","seeds")

    #set seed properties
    # cmd.set_color('palette_pink',[239/255,71/255,11/255])
    # #cmd.set_color('palette_yellow',[255/255,209/255,102/255])
    # cmd.set_color('palette_green',[6/255,214/255,160/255])
    # cmd.set_color('palette_lblue',[17/255,138/255,178/255])
    # cmd.set_color('palette_dblue',[7/255,59/255,76/255])
    # colors = ['palette_pink','paleyellow','palette_green','palette_lblue','palette_dblue']
    # cmd.set_color('palette2_green',[x/255 for x in [141,211,199]])
    # cmd.set_color('palette2_yellow',[x/255 for x in [255,255,179]])
    # cmd.set_color('palette2_purple',[x/255 for x in [190,186,218]])
    # cmd.set_color('palette2_red',[x/255 for x in [251,128,114]])
    # cmd.set_color('palette2_blue',[x/255 for x in [128,177,211]])
    # colors = ['palette2_green','palette2_yellow','palette2_purple','palette2_red','palette2_blue']
    cmd.set_color('palette3_red',[x/255 for x in [228,26,28]])
    cmd.set_color('palette3_blue',[x/255 for x in [55,126,184]])
    cmd.set_color('palette3_green',[x/255 for x in [77,175,74]])
    cmd.set_color('palette3_purple',[x/255 for x in [152,78,163]])
    cmd.set_color('palette3_orange',[x/255 for x in [250,175,0]])
    colors = ['palette3_red','palette3_blue','palette3_green','palette3_orange','palette3_purple']

    cmd.set("stick_transparency",0.25,"seeds")
    cmd.set("stick_transparency",0,"seeds")
    cmd.set("stick_radius",0.07,"seeds")

    group_seeds = cmd.get_object_list('fused_seeds')
    # group_seeds = cmd.get_object_list('seeds')
    for i,obj in enumerate(group_seeds):
        cmd.color(colors[i%len(colors)],obj)

    #set peptide properties
    cmd.color("gray80","peptide")
    cmd.set("stick_radius",0.21,"peptide")
    cmd.set("stick_transparency",0.45,"peptide")

    # cmd.set("transparency_mode",2)

    return
