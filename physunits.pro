pro physunits
  u={unit,value:0.d0,unit:'',name:''}

  PhU={!physunits, $
;Physical constants
       C    :{unit,2.99792458d10, 'cm/s',             'Speed of light'},         $
       G    :{unit,6.67259d-8,    'dyn*cm^2/g^2',     'Gravitational constant'}, $
       h    :{unit,6.6260755d-27, 'erg*s',            'Planck''s constant'}, $
       e    :{unit,4.8032068d-10, 'esu',              'Electron charge'}, $
       m_e  :{unit,9.1093897d-28, 'g',                'Mass of electron'}, $
       m_p  :{unit,1.6726231d-24, 'g',                'Mass of proton'}, $
       amu  :{unit,1.6605402d-24, 'g',                'Atomic mass unit'}, $
       R    :{unit,109737.31534d0,'1/cm',             'Rydberg constant'}, $
       k    :{unit,1.380658d-16,  'erg/K',            'Boltzmann constant'}, $
       sigma:{unit,5.670400d-5,   'erg/(cm^2*s*K^4)', 'Stefan-Boltzmann constant'}, $
;Astronomical constants
       R_earth:{unit,6378.14d5,   'cm',    'Earth radius'}, $
       M_earth:{unit,5.9737d27,   'g',     'Earth mass'}, $
       R_sun  :{unit,6.9599d10,   'cm',    'Solar radius'}, $
       L_sun  :{unit,3.826d33,    'erg/s', 'Solar luminosity'}, $
       M_sun  :{unit,1.989d33,    'g',     'Solar mass'}, $
       AU     :{unit,1.4959787d13,'cm',    'Astronomical unit'}, $
       Year   :{unit,3.1556909d7, 's',     'Year'}, $
       LY     :{unit,9.4605233d17,'cm',    'Light year'}, $
       parsec :{unit,3.085678d18, 'cm',    'parsec'}  $
      }
  defsysv,'!PhU',PhU
end
