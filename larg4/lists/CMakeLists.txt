include_directories($ENV{GEANT4_FQ_DIR}/include)
include_directories($ENV{XERCES_C_INC})

#include_directories(./include)
#include_directories(models)
#include_directories(builders)
#include_directories(processes)

#set(ENV{G4CASCADE_DEBUG_INTERFACE} 1)
add_definitions( -DG4CASCADE_DEBUG_INTERFACE=1 )
add_definitions( -DMYG4CASCADE_DEBUG_INTERFACE=1 )

art_make(
  LIB_LIBRARIES
    nug4_G4Base
    nug4_ParticleNavigation
    ${ART_UTILITIES}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    canvas
    cetlib_except
    ${ROOT_CORE}
    ${ROOT_HIST}
    ${ROOT_PHYSICS}
    ${ROOT_GEOM}
    ${ROOT_MATHCORE}
    ${CLHEP}
    ${G4DIGITS_HITS}
    ${G4EVENT}
    ${G4GEOMETRY}
    ${G4GLOBAL}
    ${G4GRAPHICS_REPS}
    ${G4MATERIALS}
    ${G4INTERCOMS}
    ${G4PARTICLES}
    ${G4PHYSICSLISTS}
    ${G4PROCESSES}
    ${G4TRACK}
    ${G4RUN}
  SERVICE_LIBRARIES
    ${ART_FRAMEWORK_CORE}
    ${MF_MESSAGELOGGER}
)

install_headers()
install_source()

# install macro files in lib
#file(GLOB g4_macros *.mac )
#install_fw( LIST ${g4_macros} )

