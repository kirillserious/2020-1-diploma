function h=myfigure(width,height)
 % создаёт figure с заданной шириной и высотой в см
 
 if nargin<2 % если дана только ширина, использует золотое сечение для высоты
     height = width * 2/(1+sqrt(5)); % высота=ширина/золотое сечение
 end
 
 set(0,'units','centimeters')
 scrsz=get(0,'screensize'); % размер экрана в см
 % положение и размер картинки
 position=[scrsz(3)/2-width/2 scrsz(4)/2-height/2 width height];
 h=figure;
 set(h,'units','centimeters')
 % устанавливаем размер
 set(h,'position',position)
 % Set screen and figure units back to pixels
 set(0,'units','pixel')
 set(h,'units','pixel')

end