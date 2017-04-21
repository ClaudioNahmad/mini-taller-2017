
Base= "../Gadget2/output/"  ;;;; path where stuff is stored

Num=  3                ;;;;; number of snapshot


exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strcompress(strmid(exts,strlen(exts)-3,3),/remove_all)

f=base+"groups_catalogue/fof_special_catalogue_"+exts     ;;;; read in group catalogue
print, f
openr,1,f
Ngroups=0L
readu,1,Ngroups
print,"Ngroups= ", Ngroups
GroupLen=lonarr(Ngroups)                            ;;;; Table with group lengths
GroupOffset=lonarr(Ngroups)
GroupMass= fltarr(Ngroups) 
GroupCM= fltarr(3, Ngroups)
GroupNspecies= lonarr(3, Ngroups)
GroupMspecies= fltarr(3, Ngroups)
GroupSfr= fltarr(Ngroups)
GroupMetallicities= fltarr(2, Ngroups)
Mcold= fltarr(Ngroups) 
SigmaStars= fltarr(Ngroups) 
SigmaDM= fltarr(Ngroups) 
readu,1, GroupLen, GroupOffset, GroupMass, GroupCM, GroupNspecies, $
         GroupMspecies, GroupSfr, GroupMetallicities
readu,1, Mcold
;readu,1, SigmaStars, SigmaDM
close,1
 

Ngas=   GroupNSpecies(0,*)
Ndm=    GroupNSpecies(1,*)
Nstars= GroupNspecies(2,*)

Mgas=   GroupMSpecies(0,*)
Mdm=    GroupMSpecies(1,*)
Mstars= GroupMspecies(2,*)

Zgas=   GroupMetallicities(0,*)
Zstars= GroupMetallicities(1,*)

print
print, "Largest group has ", GroupLen(0), " particles."
print, "Gas:      ", Ngas(0)
print, "Dm:       ", Ndm(0)
print, "Stars:    ", Nstars(0)
print, "<Zgas>=   ", Zgas(0)
print, "<Zstars>= ", Zstars(0)
print

end
