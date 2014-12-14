#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "parsecfg.h"
#include "common.h"

typedef struct ConfigEntry
{
    char *name; // имя параметра
    char *value; // значение параматера
} ConfigEntry;

struct ConfigFile
{
    ConfigEntry *v; // массив параметров
    int used; // количество элементов в массиве параметров
    int alc; // количество блоков, выделенных под массив параметров
};

int 
sort_func(const void *p1, const void *p2) {
    ConfigEntry *conf1 = (ConfigEntry *) p1;
    ConfigEntry *conf2 = (ConfigEntry *) p2;
    return strcmp(conf1->name, conf2->name);
}

ConfigFile *config_free(ConfigFile *config)
{
    if (config != NULL) {
        int i;
        for (i = 0; i < config->used; i++) {
            free(config->v[i].name);
            free(config->v[i].value);
        }
        free(config->v);
        config->v = NULL;
        free(config);
        config = NULL;
    }
    return NULL;
}

void
error_invalid_chr(const char *file, int line, char chr)
{
    fprintf(stderr, "config_read: %s: on line %d - invalid \
    character '%c'\n", file, line, chr);
}

static void
parse_error(const char *file, int line)
{
    fprintf(stderr, "Syntax error in line %d of %s\n", line, file);
}

ConfigFile *
config_read(const char *path)
{
    FILE *in = NULL;
    ConfigFile *cfg = NULL;
    char *k = NULL, *v = NULL, *str = NULL;
    if (!(in = fopen(path, "r"))) {
        error_open(path);
        goto fail;
    }
    cfg = (ConfigFile *) calloc(1, sizeof(*cfg));
    cfg->used = 0;
    cfg->alc = 1;
    cfg->v = (ConfigEntry *) calloc(1, sizeof(cfg->v[0]));
    cfg->used = 0;
    int line = 0, len = 0;
    while ((str = getline2(in)) != NULL) {
        line++;
        char *p = NULL;
        if ((p = strchr(str, '#'))) {
            *p = 0;
        }
        p = str;
        len = strlen(p);
        for (; isspace(p[len - 1]) && len > 0; len--) {}
        if (len <= 0) {
            free(str);
            continue;
        }
        p[len] = '\0';
        k = str;
        for (; isspace(*k); k++) {}
        if (*k == '\0' || *k == '\n') {
            free(str);
            continue;
        }
        v = k;
        if (!isalpha(*v) && *v != '_') {
            parse_error(path, line);
            goto fail;
        }
        for (; isalpha(*v) || isdigit(*v) || 
               *v == '_' || *v == '-'; v++){}
        if (!isspace(*v) && *v != '=') {
            parse_error(path, line);
            goto fail;
        }
        if (*v != '=') {
            *v = '\0';
            v++;
            for (; isspace(*v); v++){}
            if (*v != '=') {
                parse_error(path, line);
                goto fail;
            }
        } else {
            *v = '\0';
        }
        v++;
        for (; isspace(*v); v++){}
        char *cur = v;
        for (; *cur != '\n' && *cur != '\0'; cur++){}
        *cur = '\0';
        if (cfg->used == cfg->alc) {
            cfg->alc *= 2;
            cfg->v = (ConfigEntry *) realloc(cfg->v, 
                cfg->alc * sizeof(cfg->v[0]));
        }
        cfg->v[cfg->used].name = strdup(k);
        cfg->v[cfg->used].value = strdup(v);
        free(str);
        str = NULL;
        cfg->used++;
    }
    fclose(in);
    in = NULL;
    qsort(cfg->v, cfg->used, sizeof(cfg->v[0]), sort_func);
    int i;
    for (i = 1; i < cfg->used; i++) {
        if (!strcmp(cfg->v[i].name, cfg->v[i - 1].name)) {
            fprintf(stderr, "Duplicate parameter %s in %s\n", 
                cfg->v[i].name, path);
            goto fail;
        }
    }
    
    return cfg;
fail:
    if (in) {
        fclose(in);
        in = NULL;
    }
    cfg = config_free(cfg);
    free(str);
    return NULL;
}

const char *
config_get(ConfigFile *config, const char *name)
{
    if (!name || !config) {
        return NULL;
    }
    int left = 0, right = config->used - 1, med;
    while (left <= right) {
        med = (right + left) / 2;
        int fl = strcmp(config->v[med].name, name);
        if (fl < 0) {
            left = med + 1;
        } else if (fl > 0) {
            right = med - 1;
        } else {
            return config->v[med].value;
        }
    }
    return NULL;
}

int 
config_get_int(
    ConfigFile *config,
    const char *name,
    int *p_int)
{
    if (!name || !config || !p_int) {
        return 0;
    }
    const char *value = config_get(config, name);
    if (!value) {
        return 0;
    }
    char *endp = NULL;
    errno = 0;
    int num = strtol(value, &endp, 10);
    if (errno || *endp) {
        return -1;
    }
    *p_int = num;
    return 1;
}

void
config_print(ConfigFile *config, FILE *f)
{
    if (!config) {
        return;
    }
    int i;
    if (!f) {
        f = stdout;
    }
    for (i = 0; i < config->used; i++) {
        fprintf(f, "%s = \"%s\"\n", config->v[i].name, config->v[i].value);
    }
}
