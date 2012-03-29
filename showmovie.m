function showmovie(cinemri,ps)

  [n1,n2,n3]=size(cinemri);

  mincine=min(cinemri(:));
  maxcine=max(cinemri(:));

  for k=1:n3
    clf;
    imagesc(cinemri(:,:,k));
    colormap(gray);
    caxis([mincine maxcine]);
    axis('image')
    colorbar;
    title(num2str(k));
    pause(ps);
  end

