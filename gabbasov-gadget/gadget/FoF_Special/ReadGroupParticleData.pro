
Base= "../Gadget2/output/"  ;;;; path where stuff is stored

Num = 3   ;;;; number of dump

GrNr = 5   ;;; numer of group for which the particle data should be extracted


exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strcompress(strmid(exts,strlen(exts)-3,3),/remove_all)

f=base+"groups_catalogue/fof_special_catalogue_"+exts     ;;;; read in group catalogue
print, f
openr,1,f
Ngroups=0L
readu,1,Ngroups
print,"Ngroups= ", Ngroups
GroupLen=lonarr(Ngroups)                        ;;; Table with group lengths
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




f=base+"groups_particles/fof_special_particles_"+exts     ;;;; read in group catalogue
print, f
openr,1,f
Ntot=0L
readu,1,Ntot
print,"Ntot= ", Ntot
ParticlesPos = fltarr(3, Ntot)
readu,1, ParticlesPos
close,1



N = GroupLen(GrNr)

; extract the particles of group number GrNr

PP = ParticlesPos(*, GroupOffset(GrNr):GroupOffset(GrNr) + N -1)



end
