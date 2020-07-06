######################################
#
# Makefile for cabsim
#
######################################

# where to find the source code - locally in this case
CABSIM_SITE_METHOD = local
CABSIM_SITE = $($(PKG)_PKGDIR)/

# even though this is a local build, we still need a version number
# bump this number if you need to force a rebuild
CABSIM_VERSION = 1

# dependencies (list of other buildroot packages, separated by space)
CABSIM_DEPENDENCIES =

# LV2 bundles that this package generates (space separated list)
CABSIM_BUNDLES = cabsim.lv2

# call make with the current arguments and path. "$(@D)" is the build directory.
CABSIM_TARGET_MAKE = $(TARGET_MAKE_ENV) $(TARGET_CONFIGURE_OPTS) $(MAKE) -C $(@D)/source


# build command
define CABSIM_BUILD_CMDS
	$(CABSIM_TARGET_MAKE)
endef

# install command
define CABSIM_INSTALL_TARGET_CMDS
	$(CABSIM_TARGET_MAKE) install DESTDIR=$(TARGET_DIR)
endef


# import everything else from the buildroot generic package
$(eval $(generic-package))
