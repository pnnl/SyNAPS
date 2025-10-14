%% function to transform Ai into new variable space
function A_xfrm = transform_Ai(Ai,Tmap)

pvar = length(Ai);

for ivar = pvar:-1:1
    A_xfrm{ivar} = Ai{1} * Tmap(1,ivar);
    for jvar = 2:pvar
        A_xfrm{ivar} = A_xfrm{ivar} + Ai{jvar} * Tmap(jvar,ivar);
    end
end

end