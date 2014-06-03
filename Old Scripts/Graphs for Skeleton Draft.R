#####
# Figure 1 - Taxonomy Space
#####

par(mfrow = c(3, 2))
par(mar = c(5, 5, 4, 2) + 0.1)
MRM.point.plot("J", "bsor", "taxonomic")
par(usr = c(1,10,1,10))
text("a)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.point.plot("S", "bsor", "taxonomic")
par(usr = c(1,10,1,10))
text("b)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.point.plot("J", "bsim", "taxonomic")
par(usr = c(1,10,1,10))
text("c)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.point.plot("S", "bsim", "taxonomic")
par(usr = c(1,10,1,10))
text("d)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.point.plot("J", "bnes", "taxonomic")
par(usr = c(1,10,1,10))
text("e)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.point.plot("S", "bnes", "taxonomic")
par(usr = c(1,10,1,10))
text("f)", x = 1.5, y = 9.5)

#####
# Figure 2 - Function Space
#####

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.point.plot("J", "bsor", "functional")
par(usr = c(1,10,1,10))
text("a)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.point.plot("S", "bsor", "functional")
par(usr = c(1,10,1,10))
text("b)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.point.plot("J", "bsim", "functional")
par(usr = c(1,10,1,10))
text("c)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.point.plot("S", "bsim", "functional")
par(usr = c(1,10,1,10))
text("d)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.point.plot("J", "bsne", "functional")
par(usr = c(1,10,1,10))
text("e)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.point.plot("S", "bsne", "functional")
par(usr = c(1,10,1,10))
text("f)", x = 1.5, y = 9.5)


#####
# Figure 3 - Taxonomy Time
#####

par(mfrow = c(3, 2))
par(mar = c(5, 5, 4, 2) + 0.1)
MRM.time.point.plot("J", "bsor", "taxonomic")
par(usr = c(1,10,1,10))
text("a)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.time.point.plot("S", "bsor", "taxonomic")
par(usr = c(1,10,1,10))
text("b)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.time.point.plot("J", "bsim", "taxonomic")
par(usr = c(1,10,1,10))
text("c)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.time.point.plot("S", "bsim", "taxonomic")
par(usr = c(1,10,1,10))
text("d)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.time.point.plot("J", "bnes", "taxonomic")
par(usr = c(1,10,1,10))
text("e)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.time.point.plot("S", "bnes", "taxonomic")
par(usr = c(1,10,1,10))
text("f)", x = 1.5, y = 9.5)

#####
# Function Time
#####

par(mfrow = c(3, 2))
par(mar = c(5, 5, 4, 2) + 0.1)
MRM.time.point.plot("J", "bsor", "functional")
par(usr = c(1,10,1,10))
text("a)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.time.point.plot("S", "bsor", "functional")
par(usr = c(1,10,1,10))
text("b)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.time.point.plot("J", "bsim", "functional")
par(usr = c(1,10,1,10))
text("c)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.time.point.plot("S", "bsim", "functional")
par(usr = c(1,10,1,10))
text("d)", x = 1.5, y = 9.5)

par(mar = c(5, 5, 4, 2) + 0.1)
MRM.time.point.plot("J", "bsne", "functional")
par(usr = c(1,10,1,10))
text("e)", x = 1.5, y = 9.5)

par(mar = c(5, 4, 4, 2) + 0.1)
MRM.time.point.plot("S", "bsne", "functional")
par(usr = c(1,10,1,10))
text("f)", x = 1.5, y = 9.5)